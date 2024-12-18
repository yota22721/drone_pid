#server
from os import terminal_size
import socket
server_socket = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
server_socket.bind(("127.0.0.1",45621))
server_socket.listen(1)
client_socket,client_address = server_socket.accept()
print("Connection from",client_address)
while True:
  data = client_socket.recv(1024)
  str_data = data.decode()
  if not data : break
  print("Received from client",str_data)
  client_socket.send(data)
#client_address.close()
client_socket.close()