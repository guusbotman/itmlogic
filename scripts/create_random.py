import pickle
import random
import json
from pyproj import Transformer

from shapely import Point

# Set the boundaries of the Copernicus_DSM_30_N51_W001 tile
lat_min, lat_max = 51.0, 52.0
lon_min, lon_max = -1.0, 0.0


class receivers:
    def __init__(self, lat_min, lat_max, lon_min, lon_max):
        self.wgs = []
        self.json = []
        self.linear = []
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)

        for i in range(1000):
            lon = random.uniform(lon_min, lon_max)
            lat = random.uniform(lat_min, lat_max)
            x, y = transformer.transform(lon, lat)

            receiver = {
                'type': 'Feature',
                'geometry': {
                    'type': 'Point',
                    'coordinates': (
                        lon,  # Random longitude within the range
                        lat   # Random latitude within the range
                    )
                },
                'properties': {
                    'id': f'receiver_{i+1}'
                }
            }
            self.json.append(receiver)
            self.wgs.append(Point(lon, lat))
            self.linear.append(Point(x, y))

# # If you want to save this list to a file (optional)
# with open('receivers.json', 'w') as f:
#     json.dump(receivers, f, indent=4)

def main():
    receivers_1 = receivers(lat_min, lat_max, lon_min, lon_max)
    with open('receivers.pkl', 'wb') as f:
        pickle.dump(receivers_1, f)

if __name__ == "__main__":
    main()


