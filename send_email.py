# Install Courier SDK: pip install trycourier
from trycourier import Courier 
import os
import sys

#try:
	#if os.environ['EMAIL_CONTENT']:
    	#	content = os.environ['EMAIL_CONTENT']
#except KeyError:
#	content = "The Tsallis machine hsa finished!!!"
file1 = open("tmpf95.txt","r")#write by fortran
content = file1.read()
file1.close()

client = Courier(auth_token="token generated by courier account")

resp = client.send_message(
  message={
    "to": {
      "email": "email address",
    },
    "template": "YBKC9FHVGD4AY2MBYDDHP8HYS2TH",
    "data": {
      "recipientName": "Hermes",
      "content":content,
    },
  }
)

print(resp['requestId'])
