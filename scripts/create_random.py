import random
import json

# Set the boundaries of the Copernicus_DSM_30_N51_W001 tile
lat_min, lat_max = 51.0, 52.0
lon_min, lon_max = -1.0, 0.0

# Generate a list of 1000 random receivers within the boundaries
receivers = []
for i in range(1000):
    receiver = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': (
                random.uniform(lon_min, lon_max),  # Random longitude within the range
                random.uniform(lat_min, lat_max)   # Random latitude within the range
            )
        },
        'properties': {
            'id': f'receiver_{i+1}'
        }
    }
    receivers.append(receiver)

# If you want to save this list to a file (optional)
with open('receivers.json', 'w') as f:
    json.dump(receivers, f, indent=4)
