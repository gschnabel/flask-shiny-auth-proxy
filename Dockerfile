FROM rocker/shiny:latest

# install R packages required by the Shiny app
RUN Rscript -e "install.packages(c('shinydashboard', 'htmlTable', 'DT', 'devtools', 'MCMCpack', 'mvtnorm', 'ellipse', 'shiny'))"

# add app by copying it into the image
# (alternatively the directory of the app may
# be mapped into the container during container startup)
COPY ./shiny-app/* /srv/shiny-server/

# disable iframe-eventsource and iframe-htmlfile
# communication protocols to avoid that app greys out
RUN echo 'sanitize_errors off;disable_protocols iframe-eventsource iframe-htmlfile;' >> /etc/shiny-server/shiny-server.conf

EXPOSE 3838
