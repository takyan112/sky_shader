/* https://blog.uhawkvr.com/rendering/rendering-volumetric-clouds-using-signed-distance-fields/
 * Rendering volumetric clouds using signed distance fields
 *
 * https://www.shadertoy.com/view/slSXRW
 * Production Sky Rendering --- Created by AndrewHelmer
 */

shader_type sky;
render_mode use_debanding;


/**
 * 128x32x128 Worley noise
 */
uniform sampler3D cloudNoise;

/**
 * transmittance lookup table
 */
uniform sampler2D transmittanceLUT;

/**
 * Controls the overall density of clouds in the sky.
 * 0.0 = clear sky, 1.0 = fully overcast.
 */
uniform float cloudCoverage: hint_range(0.0, 1.0) = 0.5;

// Units are in megameters.
const float groundRadiusMM = 6.360;
const float atmosphereRadiusMM = groundRadiusMM + 0.1; // 100km

const float cloudTopMM = groundRadiusMM + 0.0120; // 12km
const float cloudBottomMM = groundRadiusMM + 0.003; // 3km

// These are per megameter.
const vec3 rayleighScatteringBase = vec3(5.802, 13.558, 33.1);
const vec3 rayleighAbsorptionBase = vec3(0.0);
const vec3 rayleighExtinctionBase = rayleighScatteringBase + rayleighAbsorptionBase;

const float mieScatteringBase = 3.996;
const float mieAbsorptionBase = 4.4;
const float mieExtinctionBase = mieScatteringBase + mieAbsorptionBase;

const vec3 ozoneScatteringBase = vec3(0.0);
const vec3 ozoneAbsorptionBase = vec3(0.650, 1.881, .085);
const vec3 ozoneExtinctionBase = ozoneScatteringBase + ozoneAbsorptionBase;

// Steps
const int raymarchScatteringSteps = 8;


// Cornette-Shanks phase function
float getMiePhase(float cosTheta) {
	const float g = 0.8;
	const float scale = 3.0/(8.0*PI);

	float num = (1.0-g*g)*(1.0+cosTheta*cosTheta);
	float denom = (2.0+g*g)*pow((1.0 + g*g - 2.0*g*cosTheta), 1.5);

	return scale*num/denom;
}


float getRayleighPhase(float cosTheta) {
	const float k = 3.0/(16.0*PI);
	return k*(1.0+cosTheta*cosTheta);
}


// https://research.nvidia.com/labs/rtr/approximate-mie/publications/approximate-mie.pdf
float getCloudPhase(float cosTheta) {
	const float d = 15.0; // water droplet diameter in µm
	const float hg_g = exp(-0.0990567/(d-1.67154));
	const float d_g = exp(-2.20679/(d+3.91029)-0.428934);
	const float a = exp(3.62489-8.29288/(d+5.52825));
	const float w_d = exp(-0.599085/(d-0.641583)-0.665888);

	float hg = 1.0/(4.0*PI)*(1.0-hg_g*hg_g)/pow(1.0+hg_g*hg_g-2.0*hg_g*cosTheta, 1.5);
	float draine = 1.0/(4.0*PI)*(1.0-d_g*d_g)/pow(1.0+d_g*d_g-2.0*d_g*cosTheta, 1.5)*(1.0+a*cosTheta*cosTheta)/(1.0+a*(1.0+2.0*d_g*d_g)/3.0);

	return (1.0 - w_d) * hg + w_d * draine;
}


void getScatteringValues(vec3 pos,
						 out vec3 rayleighScattering,
						 out float mieScattering,
						 out vec3 extinction) {
	float altitudeKM = (length(pos)-groundRadiusMM)*1000.0;
	float rayleighDensity = exp(-altitudeKM/8.0);
	float mieDensity = exp(-altitudeKM/1.2);
	float ozoneDensity = max(0.0, 1.0 - abs(altitudeKM-25.0)/15.0);

	rayleighScattering = rayleighScatteringBase*rayleighDensity;

	mieScattering = mieScatteringBase*mieDensity;

	extinction = rayleighExtinctionBase*rayleighDensity + mieExtinctionBase*mieDensity + ozoneExtinctionBase*ozoneDensity;
}


bool isIntersectSphere(vec3 ro, vec3 rd, float rad) {
	float b = dot(ro, rd);
	float c = dot(ro, ro) - rad*rad;
	return (c <= 0.0 || b <= 0.0) && (b*b >= c);
}

// return (distance to sphere, distance through sphere)
vec2 rayIntersectSphere(vec3 ro, vec3 rd, float rad) {
	float b = dot(ro, rd);
	float c = dot(ro, ro) - rad*rad;
	if (c > 0.0 && b > 0.0) return vec2(1e32, 0.0);
	float d = b*b - c;
	if (d < 0.0) return vec2(1e32, 0.0);
	float s = sqrt(d);
	if (c < 0.0) return vec2(0.0, -b + s);
	return vec2(-b - s, s + s);
}


vec3 getSunTransmittance(vec3 pos, vec3 sunDir) {
	float height = length(pos);
	float height_relative = (height - groundRadiusMM)/(atmosphereRadiusMM - groundRadiusMM);
	vec3 up = pos / height;
	float cos_zenith = dot(sunDir, up);
	vec2 uv = vec2(clamp(0.5 + 0.5*cos_zenith, 0.0, 1.0), clamp(height_relative, 0.0, 1.0));
	return texture(transmittanceLUT, uv).rgb;
}


vec3 raymarchScattering(vec3 pos, vec3 rayDir, vec3 sunDir) {
	vec2 hitInfoAtmosphere = rayIntersectSphere(pos, rayDir, atmosphereRadiusMM);
	vec2 hitInfoGround = rayIntersectSphere(pos, rayDir, groundRadiusMM);

	float dstToAtmosphere = hitInfoAtmosphere.x;
	float dstThroughAtmosphere = min(hitInfoAtmosphere.y, hitInfoGround.x - dstToAtmosphere);

	if (dstThroughAtmosphere <= 0.0) return vec3(0.0);

	float stepSize = dstThroughAtmosphere / float(raymarchScatteringSteps);

	float cosTheta = dot(rayDir, sunDir);

	float miePhaseValue = getMiePhase(cosTheta);
	float rayleighPhaseValue = getRayleighPhase(cosTheta);

	vec3 lum = vec3(0.0);
	vec3 transmittance = vec3(1.0);
	pos += dstToAtmosphere*rayDir;

	for (int i = 0; i < raymarchScatteringSteps; i++) {
		float dt = (i == 0 ? 0.3 : 1.0) * stepSize;
		pos += dt*rayDir;

		vec3 rayleighScattering, extinction;
		float mieScattering;
		getScatteringValues(pos, rayleighScattering, mieScattering, extinction);

		vec3 sampleTransmittance = exp(-dt*extinction);

		vec3 sunTransmittance = getSunTransmittance(pos, sunDir);

		const vec3 psiMS = vec3(0.0); // multiple-scattering here

		vec3 rayleighInScattering = rayleighScattering*(rayleighPhaseValue*sunTransmittance + psiMS);
		vec3 mieInScattering = mieScattering*(miePhaseValue*sunTransmittance + psiMS);
		vec3 inScattering = (rayleighInScattering + mieInScattering);

		// Integrated scattering within path segment.
		vec3 scatteringIntegral = (inScattering - inScattering * sampleTransmittance) / extinction;

		lum += scatteringIntegral*transmittance;

		transmittance *= sampleTransmittance;
	}
	return lum;
}


vec3 jodieReinhardTonemap(vec3 c){
	// From: https://www.shadertoy.com/view/tdSXzD
	float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
	vec3 tc = c / (c + 1.0);
	return mix(c / (l + 1.0), tc, tc);
}


vec3 sunWithBloom(vec3 rayDir, vec3 sunDir) {
	const float sunSolidAngle = 0.53*PI/180.0;
	const float minSunCosTheta = cos(sunSolidAngle);

	float cosTheta = dot(rayDir, sunDir);
	if (cosTheta >= minSunCosTheta) return vec3(1.0);

	float offset = minSunCosTheta - cosTheta;
	float gaussianBloom = exp(-offset*50000.0)*0.5;
	float invBloom = 1.0/(0.02 + offset*300.0)*0.01;
	return vec3(gaussianBloom+invBloom);
}


float sdf_smoothIntersection(float d1, float d2, float k)
{
	k *= 4.0;
	float h = max(k-abs(d1-d2),0.0);
	return max(d1, d2) + h*h*0.25/k;
}


float sdf_cloud(vec3 pos, float time) {
	float len = length(pos);
	float shape = max(len - cloudTopMM, cloudBottomMM - len);

	float cloud = texture(cloudNoise, pos * vec3(5.0, 20.0, 5.0) + time * vec3(0.0003, 0.0, 0.0002)).r;
	float noise = texture(cloudNoise, pos * vec3(100.0, 400.0, 100.0) + time * vec3(0.01, 0.0, 0.03)).r;

	return sdf_smoothIntersection(shape, (cloud - cloudCoverage) * 0.01 - (noise - cloudCoverage) * 0.0005, 0.001);
}


float cloud_density(float sdf) {
	return clamp(sdf * -1e7, 0.0, 200.0);
}


float getCloudTransmittance(vec3 pos, vec3 sunDir, float time) {
	if (isIntersectSphere(pos, sunDir, groundRadiusMM)) {
		return 0.0;
	}
	float transmittance = 1.0;
	const float dt = 0.002; // 2000m
	for (int i = 0; i < 4; i++) {
		float sdf = sdf_cloud(pos, time);
		if (sdf < 0.0) {
			float density = cloud_density(sdf);
			transmittance *= exp(-dt*density);
		}
		pos += (sdf < 0.0 ? dt : max(sdf, dt)) * sunDir;
	}
	return transmittance;
}

vec4 raymarchCloud(vec3 pos, vec3 rayDir, vec3 sunDir, float time) {
	bool hitTop = isIntersectSphere(pos, rayDir, cloudTopMM);
	bool hitGnd = isIntersectSphere(pos, rayDir, groundRadiusMM);
	float len = length(pos);
	if (len < groundRadiusMM ||
		(len > cloudTopMM && !hitTop) ||
		(len < cloudBottomMM && hitGnd)) {
		return vec4(0.0, 0.0, 0.0, 1.0);
	}
	if (len < cloudBottomMM) {
		pos += rayIntersectSphere(pos, rayDir, cloudBottomMM).x * rayDir;
	} else if (len > cloudTopMM) {
		pos += rayIntersectSphere(pos, rayDir, cloudTopMM).x * rayDir;
	}

	float cosTheta = dot(rayDir, sunDir);

	float miePhaseValue = getCloudPhase(cosTheta);

	vec3 lum = vec3(0.0);
	float transmittance = 1.0;
	const float dt = 0.002; // 2000m

	for (int i = 0; i < 64; i++) {
		float sdf = sdf_cloud(pos, time);
		if (sdf > 0.01) {
			break;
		}
		if (sdf < 0.0) {

			vec3 sunTransmittance = getSunTransmittance(pos, sunDir);

			sunTransmittance *= getCloudTransmittance(pos, sunDir, time);
			float density = cloud_density(sdf);

			float extinction = density;
			float scattering = density;

			float sampleTransmittance = exp(-dt*extinction);

			vec3 inScattering = scattering*(miePhaseValue*sunTransmittance);

			vec3 scatteringIntegral = (inScattering - inScattering * sampleTransmittance) / extinction;

			lum += scatteringIntegral*transmittance;

			transmittance *= sampleTransmittance;

			if (transmittance < 0.01) {
				transmittance = 0.0;
				break;
			}
		}
		pos += (sdf < 0.0 ? dt : max(sdf, dt)) * rayDir;
	}

	return vec4(lum, transmittance);
}


void sky()
{
	vec3 sunDir = normalize(LIGHT0_DIRECTION);
	vec3 rayDir = normalize(EYEDIR);

	// 200M above the ground.
	vec3 viewPos = vec3(0.0, groundRadiusMM + 0.0002, 0.0) + POSITION * 1e-6;

	vec3 lum = raymarchScattering(viewPos, rayDir, sunDir);

	// Bloom should be added at the end, but this is subtle and works well.
	vec3 sunLum = sunWithBloom(rayDir, sunDir);
	// Use smoothstep to limit the effect, so it drops off to actual zero.
	sunLum = smoothstep(0.002, 1.0, sunLum);
	if (length(sunLum) > 0.0) {
		if (isIntersectSphere(viewPos, rayDir, groundRadiusMM)) {
			sunLum *= 0.0;
		} else {
			// If the sun value is applied to this pixel, we need to calculate the transmittance to obscure it.
			sunLum *= getSunTransmittance(viewPos, sunDir);
		}
	}
	lum += sunLum;

	// Cloud
	vec4 cloud = raymarchCloud(viewPos, rayDir, sunDir, TIME);
	lum = mix(cloud.rgb, lum, cloud.a);

	// Tonemapping and gamma. Super ad-hoc, probably a better way to do this.
	lum *= 20.0;
	lum = pow(lum, vec3(1.3));
	lum /= (smoothstep(0.0, 0.2, clamp(sunDir.y, 0.0, 1.0))*2.0 + 0.15);

	lum = jodieReinhardTonemap(lum);

	lum = pow(lum, vec3(1.0/2.2));

	COLOR = lum;
}
