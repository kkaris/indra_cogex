"""Test how the API docs look like for FastAPI when using path parameters."""
from collections import Counter

from pydantic import create_model, Field
from fastapi import FastAPI
from typing import Dict, List, Type, Any, Callable, Iterable

from indra.statements import Evidence, Statement, Agent
from indra_cogex.client import queries
from inspect import signature, isfunction
from docstring_parser import parse
from functools import wraps

from indra_cogex.representation import Node


app = FastAPI()


def get_web_return_annotation(func: Callable) -> Type:
    """Get and translate the return annotation of a function."""
    # Get the return annotation
    sig = signature(func)
    return_annotation = sig.return_annotation
    if return_annotation is sig.empty:
        raise ValueError(f"Forgot to type annotate return of {func.__name__}")

    # Translate the return annotation:
    # Iterable[Node] -> List[Dict[str, Any]]
    # bool -> {func_name: bool}
    # Dict[str, List[Evidence]] -> Dict[str, List[Dict[str, Any]]] # str>int later
    # Iterable[Evidence] -> List[Dict[str, Any]]
    # Iterable[Statement] -> List[Dict[str, Any]]
    # Counter -> Dict[str, int]
    # Iterable[Agent] -> List[Dict[str, Any]]

    if return_annotation is Iterable[Node]:
        return List[Dict[str, Any]]
    elif return_annotation is bool:
        return Dict[str, bool]
    elif return_annotation is Dict[int, List[Evidence]]:
        return Dict[str, List[Dict[str, Any]]]
    elif return_annotation is Iterable[Evidence]:
        return List[Dict[str, Any]]
    elif return_annotation is Iterable[Statement]:
        return List[Dict[str, Any]]
    elif return_annotation is Counter:
        return Dict[str, int]
    elif return_annotation is Iterable[Agent]:
        return List[Dict[str, Any]]
    else:
        raise ValueError(f"Unrecognized return annotation {return_annotation}")


def fix_signature(*endpoint_params, endpoint_func: Callable):
    """Fix the signature of a function to match the API docs.

    The passed function will have the signature *args, **kwargs, we need to
    update it to its correct signature

    Inspired by:
    https://stackoverflow.com/questions/1409295/set-function-signature-in-python
    """
    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            # We hace either one or two parameters
            if len(endpoint_params) == 1:
                return f(endpoint_params[0], *args, **kwargs)
            elif len(endpoint_params) == 2:
                return f(endpoint_params[0], endpoint_params[1], *args, **kwargs)
            else:
                raise ValueError(f"{endpoint_func.__name__} has too many parameters")

        # Override signature and set the type annotations of the function arguments
        sig = signature(f)
        sig.return_annotation = get_web_return_annotation(endpoint_func)
        for param in endpoint_params:
            sig.parameters[param].annotation = str
        f.__signature__ = sig
        return wrapper


# def main():
for func_name in queries.__all__:
    if not isfunction(getattr(queries, func_name)) or func_name == "get_schema_graph":
        continue

    func = getattr(queries, func_name)
    func_sig = signature(func)
    client_param = func_sig.parameters.get("client")
    if client_param is None:
        print(f"{func_name} does not have a client parameter")
        continue

    # Get the docstring
    func_doc = parse(func.__doc__)
    func_doc_params = {doc_param.arg_name: doc_param for doc_param in func_doc.params}
    model_kwargs = {}
    for param_name, param in func_sig.parameters.items():
        if param_name == "client":
            continue
        if param_name not in func_doc_params:
            raise ValueError(f"Forgot to document {param_name} in {func_name}")
        if param.annotation is param.empty:
            raise ValueError(f"Forgot to type annotate {param_name} in {func_name}")
        # Special cases:
        # - evidence_map: Optional[List[Evidence]] -> skip
        if param_name == "evidence_map":
            continue

        documentation_param = func_doc_params[param_name]
        model_kwargs[param_name] = (
            param.annotation,
            Field(..., description=documentation_param.description),
        )
    # Get the return type
    if func_sig.return_annotation is func_sig.empty:
        raise ValueError(f"Forgot to type annotate return of {func_name}")
    if not func_doc.returns.description:
        raise ValueError(f"Forgot to document return of {func_name}")

    # Create the model
    camel = func_name.replace("_", " ").title().replace(" ", "")
    param_model = create_model(
        f"{camel}Model",
        **model_kwargs,
    )
    response_model = create_model(
        f"{camel}Response",
        **{
            func_name: (
                get_web_return_annotation(func),
                Field(..., description=func_doc.returns.description),
            )
        },
    )

    # Todo:
    #  - Create query json signature: this is tricky because the kwargs are
    #    dynamic based on the what the function we're trying to replicate
    #    looks like, i.e. if the
    #  - Create a function that care of making returns json serializable e.g.:
    #    - for bool returns: {func_name: bool}
    #    - for iterable returns: [x.to_json() for x in func_return]
    #    Try something like this:
    #    https://stackoverflow.com/a/33112180/10478812

    # Create the route function
    def route_func(*args, **kwargs):
        params = param_model(*args, **kwargs)
        return {func_name: func(**params.dict())}

    # Copy the docstring
    route_func.__doc__ = func.__doc__

    # Add the route
    # https://github.com/tiangolo/fastapi/issues/3029
    app.add_api_route(
        path=f"/{func_name}",
        endpoint=route_func,
        response_model=response_model,
        methods=["POST"],
    )


# if __name__ == '__main__':
#     main()
