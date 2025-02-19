ARG IMAGE=intersystemsdc/iris-community:preview
FROM $IMAGE

# use the root user to install packages
USER root   

# create a directory for the application     
WORKDIR /irisdev/app
RUN chown ${ISC_PACKAGE_MGRUSER}:${ISC_PACKAGE_IRISGROUP} /irisdev/app
USER ${ISC_PACKAGE_MGRUSER}

# Copy the source code
COPY . .

# install required packages
RUN pip3 install -r requirements.txt

# download the embeddings
RUN huggingface-cli download seyonec/PubChem10M_SMILES_BPE_450k --local-dir /usr/irissys/hfCache

# environment variables for embedded python
ENV LD_LIBRARY_PATH=${ISC_PACKAGE_INSTALLDIR}/bin
ENV IRISUSERNAME "SuperUser"
ENV IRISPASSWORD "SYS"
ENV IRISNAMESPACE "IRISAPP"
