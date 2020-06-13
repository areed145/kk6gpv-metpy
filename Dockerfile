FROM continuumio/miniconda3

LABEL maintainer="areed145@gmail.com"

WORKDIR /metpy-tasks
COPY . /metpy-tasks

RUN conda install cartopy
RUN conda install -c conda-forge matplotlib siphon metpy
RUN pip install pymongo==3.9.0 dnspython==1.16.0

ENV MONGODB_CLIENT 'mongodb+srv://kk6gpv:kk6gpv@cluster0-kglzh.azure.mongodb.net/test?retryWrites=true&w=majority'
ENV DISPLAY ''

CMD ["python", "metpy-tasks.py"]
