#pragma once
#ifndef _LATLON_UTM_CONVERTER_HPP_
#define _LATLON_UTM_CONVERTER_HPP_

#include <math.h>


#define RAD2DEG(rad) ((rad) * 180.0 / M_PI)
#define DEG2RAD(deg) ((deg) * M_PI / 180.0)


class LatLonUTMConverter
{
public:
	LatLonUTMConverter();
	~LatLonUTMConverter();
	LatLonUTMConverter(const LatLonUTMConverter&) = default;

public:
	// return UTM zone letter for the given latitude, return -1 if latitude is outside the UTM limits of 84N to 80S
	static char lat_to_zone_letter(double lat);
	// return UTM zone number for the given latitude/longitude, return -1 if latitude is outside the UTM limits of 180W to 180E
	static int lon_to_zone_number(double lat, double lon);

	static void LatLon_2_XY(const double lat, const double lon, double& x, double& y, char zone[4]);
	static void XY_2_LatLon(const double x, const double y, const char zone[4], double& lat, double& lon);

private:
	// earth ellipsoid model constant
	constexpr static double R = 6378137.0, E2 = 0.00669438; // Equator Radius, Squared Eccentricity
	constexpr static double K0 = 0.9996; // scale of the central meridian

	// auxiliary variables
	constexpr static char ZONE_LETTERS[] = "CDEFGHJKLMNPQRSTUVWXX";
	constexpr static double Ep2 = E2 / (1.0 - E2);
	constexpr static double E4 = E2 * E2;
	constexpr static double E6 = E2 * E4;
	constexpr static double Ec = 0.996647189330307; // sqrt(1 - E2);
	constexpr static double Eo = (1.0 - Ec) / (1.0 + Ec);
	constexpr static double Eo2 = Eo * Eo;
	constexpr static double Eo3 = Eo * Eo2;
	constexpr static double Eo4 = Eo * Eo3;
	constexpr static double Eo5 = Eo * Eo4;
	constexpr static double M1 = 1.0 - E2 / 4.0	- 3.0 * E4 / 64.0 - 5.0 * E6 / 256.0;
	constexpr static double M2 = 3.0 * E2 / 8.0	+ 3.0 * E4 / 32.0 + 45.0 * E6 / 1024.0;
	constexpr static double M3 = 15.0 * E4 / 256.0 + 45.0 * E6 / 1024.0;
	constexpr static double M4 = 35.0 * E6 / 3072.0;
	// constexpr static double P2 = 3.0 * Eo / 2.0 - 27.0 * Eo3 / 32.0;
	constexpr static double P2 = 3.0 * Eo / 2.0 - 27.0 * Eo3 / 32.0 + 269.0 * Eo5 / 512.0;
	constexpr static double P4 = 21.0 * Eo2 / 16.0 - 55.0 * Eo4 / 32.0;
	constexpr static double P6 = 151.0 * Eo3 / 96.0;
	// constexpr static double P8 = 1097.0 * Eo4 / 512.0;
	constexpr static double P8 = 1097.0 * Eo4 / 512.0 - 417.0 * Eo5 / 128.0;
};


#endif
