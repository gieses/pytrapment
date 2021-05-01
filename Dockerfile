FROM python:3.8.3-slim

# Setup a spot for the code
WORKDIR /pytrapment

# Install Python dependencies
COPY Pipfile Pipfile
RUN pip install pipenv
RUN pipenv install Pipfile --dev

COPY pytrapment pytrapment/
COPY tests tests/

CMD ["/bin/bash"]