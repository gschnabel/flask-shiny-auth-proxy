import flask
import flask_login
import requests

app = flask.Flask(__name__)
app.secret_key = 'think twice!'


login_manager = flask_login.LoginManager()
login_manager.init_app(app)

users = {'mrano@ano.co': {'password': 'supersecret'}}


class User(flask_login.UserMixin):
    pass


@login_manager.user_loader
def user_loader(email):
    if email not in users:
        return
    user = User()
    user.id = email
    return user


@login_manager.request_loader
def request_loader(request):
    email = request.form.get('email')
    if email not in users:
        return
    user = User()
    user.id = email

    user.is_authenticated = request.form['password'] == users[email]['password']
    return user


@app.route('/login', methods=['GET', 'POST'])
def login():
    if flask.request.method == 'GET':
        return '''
		   <form action='login' method='POST'>
			<input type='text' name='email' id='email' placeholder='email'/>
			<input type='password' name='password' id='password' placeholder='password'/>
			<input type='submit' name='submit'/>
		   </form>
        '''
    email = flask.request.form['email']
    if flask.request.form['password'] == users[email]['password']:
        user = User()
        user.id = email
        flask_login.login_user(user)
        return flask.redirect(flask.url_for('protected'))

    return 'Bad login'


@app.route('/protected')
@flask_login.login_required
def protected():
    return 'Logged in as: ' + flask_login.current_user.id


@app.route('/logout')
def logout():
    flask_login.logout_user()
    return 'Logged out'


@login_manager.unauthorized_handler
def unauthorized_handler():
    return 'Unauthorized'


SITE_NAME = 'http://127.0.0.1:8765'

@app.route('/nucleardata/', defaults={'path': ''}, methods=['GET', 'POST'])
@app.route('/nucleardata/<path:path>', methods=['GET', 'POST'])
def proxy(path):
    global SITE_NAME
    print(path)
    # excluded_headers = ['content-encoding', 'content-length', 'transfer-encoding', 'connection']
    excluded_headers = ['content-encoding', 'content-length', 'transfer-encoding', 'connection']
    if flask.request.method=='GET':
        querystr = flask.request.query_string.decode() 
        target_url = f'{SITE_NAME}/{path}'
        # if querystr != '':
        #     target_url += f'?{querystr}'
        resp = requests.get(f'{target_url}', params=flask.request.args.to_dict())
        # resp = requests.get(f'{target_url}')
        headers = [(name, value) for (name, value) in  resp.raw.headers.items() if name.lower() not in excluded_headers]
        response = flask.Response(resp.content, resp.status_code, headers)
        return response
    elif flask.request.method=='POST':
        querystr = flask.request.query_string.decode() 
        target_url = f'{SITE_NAME}/{path}'
        # if querystr != '':
        #     target_url += f'?{querystr}'
        content_type = flask.request.headers.get('Content-Type')
        if content_type == 'application/json':
            resp = requests.post(f'{target_url}',json=flask.request.get_json())
        elif content_type == 'application/x-www-form-urlencoded':
            resp = requests.post(f'{target_url}',data=flask.request.form.to_dict())
        else:
            raise TypeError(f'Do not know how to deal with Content-Type {content_type}')
        headers = [(name, value) for (name, value) in resp.raw.headers.items() if name.lower() not in excluded_headers]
        response = flask.Response(resp.content, resp.status_code, headers)
        return response
    elif flask.request.method=='DELETE':
        resp = requests.delete(f'{SITE_NAME}{path}').content
        response = flask.Response(resp.content, resp.status_code, headers)
        return response


if __name__ == '__main__':
    app.run()
