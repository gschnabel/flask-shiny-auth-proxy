FROM rocker/rstudio

# install some R required stuff
RUN apt-get update -y --no-install-recommends \
    && apt-get -y install -f \
       zlib1g-dev \
       libssl-dev \
       libcurl4-openssl-dev
       
# install R packages
RUN Rscript -e "install.packages(c('shinydashboard', 'htmlTable', 'DT', 'devtools', 'MCMCpack', 'mvtnorm', 'ellipse', 'shiny'))"

# install shiny
# RUN export ADD=shiny && bash /etc/cont-init.d/add
RUN /rocker_scripts/install_shiny_server.sh

# add app to the server
COPY ./shiny-app/* /srv/shiny-server/

# update the index page
# COPY index_page/index.html /srv/shiny-server/index.html
# COPY index_page/img /srv/shiny-server/img

# try to avoid greying out of the apps
# https://stackoverflow.com/questions/44397818/shiny-apps-greyed-out-nginx-proxy-over-ssl
#RUN echo 'sanitize_errors off;disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf
# RUN echo 'sanitize_errors off;disable_protocols websocket xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf
#RUN echo 'sanitize_errors off;disable_protocols websocket xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile xdr-polling xhr-polling iframe-xhr-polling;' >> /etc/shiny-server/shiny-server.conf
# NEXT LINE IS WORKING
#RUN echo 'sanitize_errors off;disable_protocols websocket xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf
RUN echo 'sanitize_errors off;disable_protocols iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf

# websocket
# xdr-streaming
# xhr-streaming
# iframe-eventsource
# iframe-htmlfile
# xdr-polling
# xhr-polling
# iframe-xhr-polling
# jsonp-polling

EXPOSE 8787

CMD ["/init"]
