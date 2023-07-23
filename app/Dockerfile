FROM python:alpine

RUN apk update \
    && apk add --no-cache gcc musl-dev libffi-dev \
    ;

COPY ./requirements.txt .

RUN pip install --upgrade pip \
    && pip install -r requirements.txt \
    ;

RUN addgroup -g 1000 appuser && adduser -u 1000 -G appuser -h /home/appuser -D appuser

USER appuser

WORKDIR /home/appuser/app