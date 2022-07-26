FROM nginx:latest

RUN apt-get update -qq && apt-get -y install apache2-utils

ENV RAILS_ROOT /var/www/barkeeper
ARG PROJECT_DOMAIN
ARG PORT
ARG SSL_PORT
ARG PUMA_PORT
ARG RAILS_ENV

WORKDIR $RAILS_ROOT

RUN mkdir log

COPY public public/
COPY nginx.conf /tmp/docker.nginx
COPY ssl.conf* /tmp/docker.ssl
RUN envsubst '${RAILS_ROOT} ${PROJECT_DOMAIN} ${PUMA_PORT} ${PORT} ${SSL_PORT}' < /tmp/docker.nginx > /etc/nginx/conf.d/default.conf
RUN envsubst '${PROJECT_DOMAIN}' < /tmp/docker.ssl > /etc/nginx/conf.d/ssl.conf

EXPOSE ${PORT}
EXPOSE ${SSL_PORT}

CMD [ "nginx", "-g", "daemon off;" ]
