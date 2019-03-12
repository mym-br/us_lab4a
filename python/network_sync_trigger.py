#!/usr/bin/env python2
# This file is in the public domain.

# Network trigger.

import socket

HOST = '10.10.184.234'
PORT = 5544

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))
s.sendall("TRIGGER\r\n")

print("Waiting for reply message...")
reply_msg = s.recv(1024)
s.close()
print("Received: {}".format(repr(reply_msg)))
