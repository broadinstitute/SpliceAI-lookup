This folder contains the [Google Cloud Run](https://cloud.google.com/run) implementation of SpliceAI and Pangolin web service APIs used by [spliceai-lookup.broadinstitute.org](https://spliceai-lookup.broadinstitute.org)  (NOTE: Cloud Run is different from Google's Cloud Functions service).

`make` sub-commands:
  
* **build** the docker image  
* **push** the docker image to Google Artifact Registry  
* **deploy** the service to Google Cloud Run, attaching `gs://spliceai-lookup-reference-data/` as a gcsfuse volume mount
* **test** run the service locally using a `docker run` command
* **test2** run the service locally using the heavier-weight `gcloud beta code dev` command which uses kubectl

To perform any of these operations, run `make <sub-command>`.
