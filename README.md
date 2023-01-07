## Prototype of Flask reverse proxy for Shiny application

The community version of [Shiny] does not come with user
authentication. User authentication can nevertheless be
achieved by putting a reverse proxy between the user and the
Shiny server instance, which only allows successfully authenticated
users to connect to the Shiny server.

This repository demonstrates how such a reverse proxy
can be implemented using the Python `Flask` web development kit. 
As in my tests Shiny applications running behind this
reverse proxy quickly greyed out and were unresponsive,
it was necessary to make small changes to the Shiny
server configuration regarding the client-server communication.
These changes are implemented in the `Dockerfile` in this repository.

## Usage

For a fast and easy installation of the Shiny server including an
app for demonstration purposes, it is recommended to install
the free-of-charge community edition of [Docker].
If Docker is installed and you've downloaded/cloned this
repository, you can create the Docker image by running
from within the root directory of this repository:
```
docker build -t shiny-app-test .
```
You may need `sudo` or administrator priviledges to run this
and the next instruction. This instruction should complete
after a couple of minutes once all necessary components
have been downloaded from the internet.

Afterwards, you can instantiate a Docker container running
the Shiny server via 
```
docker run --rm -p 8765:3838 -d shiny-app-test
```
The first number `8765` specifies under which port the server
is reachable on your host system. It can be changed, but then
it needs also be changed in the `SITE_NAME` variable defined
at the beginning of `flask_auth_proxy.py`.

Before starting the Flask reverse proxy, make sure that the Python
packages `flask`, `flask_login` and `requests` are installed.
The reverse proxy can be launched by
```
python flask_auth_proxy.py
```
which by default exposes the reverse proxy on port 5000.

The Shiny application is available under
<http://localhost:5000/shinyapp>. As unauthorized user
it is not possible to access the app. To authenticate
yourself, go to <http://localhost:5000/login> and use
the email `mrano@ano.co` and the password `supersecret`.
After successful authentication, the app is accessible.
To log out, go to <http://localhost:5000/logout>.

## Note of caution

The only purpose of this repository is to show the basic
ingredients to set up a reverse proxy for Shiny apps
using the Flask web development kit. It is certainly
not a production ready solution due the following flaws
among potentially others:

- The password is stored in clear text in the source code,
  which should never be done. Instead, a cryptographic hash
  function, such as [SHA-256], can be used to compute
  a hash string, which is stored in place of the password
  as a minimum, but more measures should be taken, such as
  [salting]. Also, the passwords are better stored in a
  separate database than in the source code.

- The connection to the Flask instance is unencrypted and
  communication, in particular credentials, are sent in
  clear text. The client-server communication should be
  protected by using the `https` protocol and a server
  certificate.

## License

This prototype is released under the [Unlicense] license,
see the accompanying `LICENSE` file for more information.

## Helpful resources

The following resources were helpful in coming up
with a working solution for the Flask reverse proxy:

- <https://medium.com/customorchestrator/simple-reverse-proxy-server-using-flask-936087ce0afb>
- <https://pythonawesome.com/flask-user-session-management/>
- <https://stackoverflow.com/questions/44397818/shiny-apps-greyed-out-nginx-proxy-over-ssl>
- <https://github.com/hemberg-lab/scmap-shiny/blob/3d3381aa9780b258e29afd6de0b08326a6571b44/Dockerfile#L20>


[Shiny]: https://shiny.rstudio.com/
[Docker]: https://www.docker.com/
[SHA-256]: https://en.wikipedia.org/wiki/SHA-2
[salting]: https://en.wikipedia.org/wiki/Salt_(cryptography)
[Unlicense]: https://unlicense.org/
