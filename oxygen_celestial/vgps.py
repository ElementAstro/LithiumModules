# coding=utf-8

"""

Copyright(c) 2022-2023 Max Qian  <lightapt.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License version 3 as published by the Free Software Foundation.
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.
You should have received a copy of the GNU Library General Public License
along with this library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
Boston, MA 02110-1301, USA.

"""

import os, sys, signal, time, datetime, configparser, re

config_file = "location.conf"
virtualgps_dev = "/tmp/vgps"


def convert_to_sexagesimal(coord):
	"""
	Convert a string of coordinates using delimiters for minutes ('),
	seconds (") and degrees (º). It also supports colon (:).
		>>> from virtualgps import convert_to_sexagesimal
		>>> convert_to_sexagesimal(u"52:08.25:1.5\"")
		52.13791666666667
		>>> convert_to_sexagesimal(u"52:08:16.5\"")
		52.13791666666667
		>>> convert_to_sexagesimal(u"52.1:02:16.5\"")
		52.13791666666667
		>>> convert_to_sexagesimal(u"52º08'16.5\"")
		52.13791666666667
	:param coord: Coordinates in string representation
	:return: Coordinates in float representation
	"""
	elements = re.split(r'[\u00ba\':\"]', coord)

	degrees = float(elements[0])
	if (len(elements) - 1) > 0:
		# Convert minutes to degrees
		degrees += float(elements[1]) / 60
	if (len(elements) - 1) > 1:
		# Convert seconds to degrees
		degrees += float(elements[2]) / 3600
	return degrees


def nmea_checksum(sentence):
    chsum = 0
    for s in sentence:
        chsum ^= ord(s)
    return hex(chsum)[2:]

def shutdown():
	try:
		os.remove(virtualgps_dev)
	except:
		pass
	os.close(master)
	os.close(slave)
	sys.exit()

def term_handler(signum, frame):
	raise KeyboardInterrupt

# register term handler
signal.signal(signal.SIGTERM, term_handler)

if __name__ == '__main__':
	# create pseudo terminal device
	master, slave = os.openpty()
	pty = os.ttyname(slave)

	# remove leftovers before setting virtual dev
	if os.path.isfile(virtualgps_dev):
		os.remove(virtualgps_dev)

	os.symlink(pty,virtualgps_dev)

	# load location data from config
	if os.path.isfile(config_file):
		config = configparser.ConfigParser()
		config.read(config_file)
		if 'latitude' in config['default'] and 'longitude' in config['default'] and 'elevation' in config['default']:
			latitude = convert_to_sexagesimal(config['default']['latitude'])
			longitude = convert_to_sexagesimal(config['default']['longitude'])
			elevation = float(config['default']['elevation'])
		else:
			# if config wrong exit
			raise KeyboardInterrupt
	else:
		# if config does not exist exit
		raise KeyboardInterrupt

	# W or E
	if latitude > 0:
		NS = 'N'
	else:
		NS = 'S'

	# N or S
	if longitude > 0:
		WE = 'E'
	else:
		WE = 'W'

	# format for NMEA
	latitude = abs(latitude)
	longitude = abs(longitude)
	lat_deg = int(latitude)
	lon_deg = int(longitude)
	lat_min = (latitude - lat_deg) * 60
	lon_min = (longitude - lon_deg) * 60
	latitude = "%02d%07.4f" % (lat_deg, lat_min)
	longitude = "%03d%07.4f" % (lon_deg, lon_min)
	while True:
		try:
			now = datetime.datetime.utcnow()
			date_now = now.strftime("%d%m%y")
			time_now = now.strftime("%H%M%S")

			# NMEA minimal sequence:
			#$GPGGA,231531.521,5213.788,N,02100.712,E,1,12,1.0,0.0,M,0.0,M,,*6A
			#$GPGSA,A,1,,,,,,,,,,,,,1.0,1.0,1.0*30
			#$GPRMC,231531.521,A,5213.788,N,02100.712,E,,,261119,000.0,W*72

			# assemble nmea sentences
			gpgga = "GPGGA,%s,%s,%s,%s,%s,1,12,1.0,%s,M,0.0,M,," % (time_now, latitude, NS, longitude, WE, elevation)
			gpgsa = "GPGSA,A,3,,,,,,,,,,,,,1.0,1.0,1.0"
			gprmc = "GPRMC,%s,A,%s,%s,%s,%s,,,%s,000.0,W" % (time_now, latitude, NS, longitude, WE, date_now)

			# add nmea checksums
			gpgga = "$%s*%s\n" % (gpgga, nmea_checksum(gpgga))
			gpgsa = "$%s*%s\n" % (gpgsa, nmea_checksum(gpgsa))
			gprmc = "$%s*%s\n" % (gprmc, nmea_checksum(gprmc))

			os.write(master, gpgga.encode())
			os.write(master, gpgsa.encode())
			os.write(master, gprmc.encode())

			time.sleep(1)
		except KeyboardInterrupt:
			shutdown()