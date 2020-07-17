TAG1 = weisburd/spliceai-lookup-hg38:v1
TAG2 = weisburd/spliceai-lookup-hg38:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
