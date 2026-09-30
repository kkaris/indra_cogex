"""Chat page app for INDRA CoGEx"""
import logging
from pathlib import Path

import flask

from indra_cogex.apps.constants import (
    LOCAL_VUE,
    STATIC_DIR,
)

logger = logging.getLogger(__name__)

# Set corresponding endpoint in indra_cogex/apps/chat_page/app/src/components/DiscoveryApp.vue
# The prefix could interfere with the endpoint set in CloudFront/S3 (and in vue.config.js)
chat_blueprint = flask.Blueprint("chat_api", __name__, url_prefix="/api_chat")

__all__ = [
    "chat_blueprint",
]

# Serve vue app locally for testing
if LOCAL_VUE:
    from flask import send_from_directory

    if (isinstance(LOCAL_VUE, str) and not Path(LOCAL_VUE).is_dir()) or isinstance(
        LOCAL_VUE, bool
    ):
        DIST = STATIC_DIR / "vue-chat" / "dist"
    else:
        DIST = Path(LOCAL_VUE)

    logger.info(f"Serving vue app locally from {DIST}")

    @chat_blueprint.route("/vue/<path:file>")
    def serve_vue(file):
        return send_from_directory(DIST.absolute().as_posix(), file)


else:
    logger.info("Serving vue app from [not implemented]")
