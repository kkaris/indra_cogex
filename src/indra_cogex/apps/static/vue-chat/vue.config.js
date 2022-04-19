module.exports = {
    chainWebpack: config => {
        config.module
            .rule("vue")
            .use("vue-loader")
            .tap(options => {
                options.compilerOptions = {
                    ...options.compilerOptions,
                    // treat any tag that contains a dash as a custom element
                    isCustomElement: tag => tag.contains("-")
                }
                return options
            });
    }
};