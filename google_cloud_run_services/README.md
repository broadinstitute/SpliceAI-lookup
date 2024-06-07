This folder contains the [Google Cloud Run](https://cloud.google.com/run) implementation of SpliceAI and Pangolin web service APIs used by [spliceai-lookup.broadinstitute.org](https://spliceai-lookup.broadinstitute.org)  (NOTE: `Cloud Run` is different from Google's `Cloud Functions` service).

The `build_and_deploy.py` script includes the follwoing commands for building docker images, updating gencode annotations, updating the SpliceAI-lookup Google Cloud Run services, and running tests:

* **build** the docker images for the SpliceAI and Pangolin services
* **update_annotations** download Gencode annotations and reprocess them into the formats used by SpliceAI and Pangolin
* **deploy** the services to Google Cloud Run
* **test** run the service locally using a `docker run` command
* **test2** run the service locally using the heavier-weight `gcloud beta code dev` command which uses kubectl
* **run** open an interactive shell inside the container
To perform any of these operations, run `python3 build_and_deploy.py <sub-command>`.
