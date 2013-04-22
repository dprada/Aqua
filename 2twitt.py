#!/usr/bin/env python

import pbkdf2
from Crypto.Cipher import AES
import base64
import os
import getpass
import tweepy
import subprocess as sp

tt=sp.Popen("git log --pretty=format:%s -n1", stdout=sp.PIPE, shell=True).stdout.read()

if tt.startswith('@tt'):

    BLOCK_SIZE = 32
    PADDING = '{'
    pad = lambda s: s + (BLOCK_SIZE - len(s) % BLOCK_SIZE) * PADDING
    EncodeAES = lambda c, s: base64.b64encode(c.encrypt(pad(s)))
    DecodeAES = lambda c, e: c.decrypt(base64.b64decode(e)).rstrip(PADDING)

    home='/home/diego/'
    fff=open(home+'.sec/tt.key','r')
    lll=fff.readlines()
    fff.close()
    salt='\x89\x88\x8b|G\xbe6\x97'
    iv='\x1ew\x00;\x01\xf8\xf0\xd1y\xba/"\xdb\xc1\xd0\x15'
    errormes='Error! The key phrase is not valid'

    pasa=False
    contador=0
    while not pasa:
        frase=getpass.getpass('Enter Key Phrase: ')
        key = pbkdf2.PBKDF2(frase, salt).read(32)
        cipher = AES.new(key, AES.MODE_CBC, iv)
        clave = EncodeAES(cipher,errormes)
        if (clave!=lll[0][:-1]):
            print errormes
            contador+=1
        else:
            pasa=True
        if contador==3:
            print '### Three attempts without a valid Key Phrase!'
            print '### Exit.'
            os.sys.exit()

    cipher = AES.new(key, AES.MODE_CBC, iv)
    ck = DecodeAES(cipher, lll[1][:-1])
    cipher = AES.new(key, AES.MODE_CBC, iv)
    cs = DecodeAES(cipher, lll[2][:-1])
    cipher = AES.new(key, AES.MODE_CBC, iv)
    at = DecodeAES(cipher, lll[3][:-1])
    cipher = AES.new(key, AES.MODE_CBC, iv)
    ats = DecodeAES(cipher, lll[4][:-1])


    #toplevel_path=sp.Popen("git rev-parse --show-toplevel", stdout=sp.PIPE, shell=True).stdout.read()
    #toplevel_dir=toplevel_path.split('/')[-1].replace('\n','')
    #branch=sp.Popen("git rev-parse --abbrev-ref HEAD", stdout=sp.PIPE, shell=True).stdout.read()
    gitkey=sp.Popen("git log --pretty=format:%H -n1", stdout=sp.PIPE, shell=True).stdout.read()
    gitabbkey=sp.Popen("git log --pretty=format:%h -n1", stdout=sp.PIPE, shell=True).stdout.read()
    #gitauthor=sp.Popen("git log --pretty=format:%an -n1", stdout=sp.PIPE, shell=True).stdout.read()
    #gitcommauth=sp.Popen("git log --pretty=format:%cn -n1", stdout=sp.PIPE, shell=True).stdout.read()
    gitsubj=sp.Popen("git log --pretty=format:%s -n1", stdout=sp.PIPE, shell=True).stdout.read()
    gitbody=sp.Popen("git log --pretty=format:%b -n1", stdout=sp.PIPE, shell=True).stdout.read()

    #hashtags="#code #"+toplevel_dir
    subj=gitsubj.replace('@tt','')
    subj=subj.lstrip()
    if subject.endswith()!='.':
        subj=subj+'.'
    if gitbody.startswith('\n'):
        gitbody=gitbody[1:]
    gitbody.replace('\n',' ')
    if len(gitbody)>0:
        if gitbody.endsith()!='.':
            gitbody=gitbody+'.'
    tweet='[git]: "'+subject+gitbody
    if len(tweet)>116:
        tweet=tweet[:114]+'... '
    webcomm=' https://github.com/dprada/Aqua/commit/'+gitkey    
    tweet+=webcomm

    auth = tweepy.OAuthHandler(ck,cs)
    auth.set_access_token(at, ats)
    api = tweepy.API(auth)

    print '--> @tt'
