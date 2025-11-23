from flask import Flask, render_template, jsonify, request
import socketio
import eventlet
import eventlet.wsgi
import json

# Initialize Flask app and Socket.IO server
app = Flask(__name__)
sio = socketio.Server(cors_allowed_origins="*")

# Store simulation data for real-time visualization
simulation_data = {}

@sio.on('connect')
def handle_connect(sid, environ):
    print('Client connected:', sid)

@sio.on('disconnect')
def handle_disconnect(sid):
    print('Client disconnected:', sid)

@sio.on('simulation_data')
def handle_simulation_data(sid, data):
    global simulation_data
    simulation_data = data
    # Broadcast data to all connected clients
    sio.emit('update_data', data, skip_sid=sid)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/data')
def get_data():
    return jsonify(simulation_data)

if __name__ == '__main__':
    # Wrap Flask app with Socket.IO middleware
    app.wsgi_app = socketio.WSGIApp(sio, app.wsgi_app)
    # Run the server
    eventlet.wsgi.server(eventlet.listen(('', 5000)), app)
