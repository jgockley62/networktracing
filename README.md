# Igraph Sandbox

### Install Git, Docker, and Docker-Compose if needed
This is configured for AWS EC-2 instance
```{bash}
sudo yum install -y git
git version

sudo amazon-linux-extras install docker
docker version

sudo curl -L https://github.com/docker/compose/releases/download/1.22.0/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
docker-compose version
```

### Setup The Git Repo
```{bash}
git clone https://github.com/jgockley62/networktracing.git
cd networktracing
```

### Start Docker
```{bash}
sudo service docker start
sudo usermod -aG docker <USR_ID>
```

### Build RStudio Images
```{bash}
#RStudio
docker image build -t network ~/networktracing/Docker/

#Cytoscape Linux VIM
docker pull biodepot/novnc-cynetworkbma

#Caddy proxy https login 
docker build -t cyto-caddy caddy/.
```
### Build the envronment object containing login credentials
This example will create the user name of both IP's as jgockley and the password as test

```{bash}
#Copy the scratch environment file into the used .env file
sudo cp ENV_Scratch .env
 
#Insert the hashed the password for the caddy image
sudo sed -i  "s/hash/$(docker run --rm -it cyto-caddy caddy hash-password -plaintext 'test')/" .env

#Insert the User name
sudo sed -i  "s/user/jgockley/" .env

#Insert the RStudio Password 
sudo sed -i  "s/pass/test/" .env

```

### Build Containers
```{bash}

docker-compose up -d

```

### Ported Browser Access
RStudio Instance Available at: https://<AWS Instance IP>:8787

Caddy Proxy Login to Access Cytoscape NoVNC VIM Available at: https://<AWS Instance IP>:6080

### Shut Containers Down
```{bash}

docker-compose down -v

```

### Deprecated code
```{bash}

docker run -v "~/networktracing/:~/networktracing/" -e USER=<USERID> -e PASSWORD=<PassWD> -d -p 8787:8787 <ImageID>

docker run -v /home/jgockley/networktracing:/home/jgockley/networktracing -d -p 6080:6080 biodepot/novnc-cynetworkbma

```

