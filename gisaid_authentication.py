""" Authenticate with GISAID credentials """

from outbreak_data import authenticate_user

def check_authentication():
    """ Authenticate with GISAID credentials """
    authenticate_user.get_authentication()

if __name__ == "__main__":
    try:
        check_authentication()
        print("Authentication already exists")
    except:
        authenticate_user.authenticate_new_user()