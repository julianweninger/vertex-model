# NOTE ------------------------------------------------------------------------
# IMPORTANT:
#     Do not forget to increment the $MODELS_BASE_IMAGE_VERSION variable in the
#     gitlab CI!
#     The version has the shape x.y, where x is incremented only for changes
#     that are not backwards-compatible
# -----------------------------------------------------------------------------

# Use Utopia base image as starting point
# NOTE: The default setting will likely not work because no tag is specified!
ARG UTOPIA_BASE_IMAGE=ccees/utopia-base
FROM ${UTOPIA_BASE_IMAGE}

LABEL maintainer="Lukas Riedel <lriedel@iup.uni-heidelberg.de>, Yunus Sevinchan <ysevinch@iup.uni-heidelberg.de>"
