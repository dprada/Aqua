import tweepy
import subprocess as sp


toplevel_path=sp.Popen("git rev-parse --show-toplevel", stdout=sp.PIPE, shell=True).stdout.read()
toplevel_dir=toplevel_path.split('/')[-1].replace('\n','')


branch=sp.Popen("git rev-parse --abbrev-ref HEAD", stdout=sp.PIPE, shell=True).stdout.read()
gitshow=sp.Popen("git log --pretty=format:%s -n1", stdout=sp.PIPE, shell=True).stdout.read()

hashtags="#code #"+toplevel_dir
tweet=hashtags+' ['+branch+']: "'+subject+'"'

consumer_key="VB5P9mEMXWeKsWiyf9g"
consumer_secret="UsM8MdYJW7mFiZwfrpIgOZZZGEM1hAX0DpGnl3bQCos"
access_token="1364106511-CjqfjLuBQWHQeBtSqQmlWBdbz3Lpjn0gSpAV7Xs"
access_token_secret="CGEJPGgRBKHZROFENOhEXheKCFYINbBY9tvm8ILaAU"

auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
auth.set_access_token(access_token, access_token_secret)

api = tweepy.API(auth)

print tweett
#api.update_status(tweett)


