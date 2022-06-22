
import socket

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server_address = ('192.168.1.206', 10000)
sock.bind(server_address)
sock.listen(1)

conn, client_address = sock.accept()

print(conn, client_address)

conn.close()
