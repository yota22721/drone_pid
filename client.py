#client
import socket

client_socket = socket.socket()
client_socket.connect(("127.0.0.1",45621))
while True:
  data = input(">")
  client_socket.send(data.encode())
  if  not data:break
  newData = client_socket.recv(1024)
  print("Received from server:" ,str(newData.decode()))
client_socket.close()