# Image for a Python 3 development environment
FROM python:3.9-slim

# Add any tools that are needed beyond Python 3.9
RUN apt-get update && \
    apt-get install -y sudo vim git zip tree curl wget jq && \
    apt-get autoremove -y && \
    apt-get clean -y

# Create a user for development
ARG USERNAME=devops
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG DOTFILEURL="https://github.com/bashgame/cuddly-journey/raw/main/dotfiles/dotfiles.tar.gz"

# Create the user with passwordless sudo privileges
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && usermod -aG sudo $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    && curl -LS $DOTFILEURL | tar -xzC /home/$USERNAME --owner=$USER_UID --group=$USER_GID

# Set up the Python development environment
WORKDIR /app
RUN python -m pip install --upgrade pip && \
    pip install --upgrade wheel

# Enable color terminal for docker exec bash
ENV TERM=xterm-256color

EXPOSE 8000

# Become a regular user for development
USER $USERNAME
