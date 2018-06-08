#!/usr/bin/env python2

# Network trigger abort.

import socket

HOST = '10.10.184.234'
PORT = 5544

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))
s.sendall("TRIGGER ABORT\r\n")
s.close()
