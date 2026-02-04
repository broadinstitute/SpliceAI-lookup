# Use a slim Python image for a smaller footprint
FROM python:3.11-slim.bookworm

# Import proxy settings from build arguments
ARG HTTPS_PROXY

# Set environment variables for proxy -> persistent in the container
ENV http_proxy=${HTTPS_PROXY}
ENV https_proxy=${HTTPS_PROXY}

# Install redis-tools if you absolutely need redis-cli inside the container
RUN apt-get update && apt-get install -y redis-tools && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements first to leverage Docker cache
COPY olis_requirements.txt .
RUN pip install --no-cache-dir -r olis_requirements.txt

# Copy the rest of the app
COPY . .

# Set environment variables (can be overridden in docker-compose)
ENV NUM_THREADS=4
ENV TIMEOUT=1800
ENV PORT=8080

# Expose the port
EXPOSE 8080

# The command to run the app
CMD gunicorn -w $NUM_THREADS -t $TIMEOUT -b 0.0.0.0:$PORT server:app
