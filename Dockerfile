FROM biomehub/encodestatsbase:0.0.1
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    openssh-client
ENV R_REMOTES_UPGRADE=never
ENV VIRTUAL_ENV_DISABLE_PROMPT=1
RUN echo "PS1='\[\e[1;38;2;231;41;138m\]${VIRTUAL_ENV:+[$(basename -- $VIRTUAL_ENV)] }\[\e[1;38;2;117;112;179m\][[\u]]\[\033[00m\]:\[\e[1;38;2;27;158;119m\]\w/\n\[\e[1;38;2;217;95;2m\]\\$\\$\[\033[00m\] '" >> ~/.bashrc
RUN R -e "devtools::install_version('targets', version = '1.3.2')"