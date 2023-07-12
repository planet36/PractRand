#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

//for use in seeding & self-tests:
#include "PractRand/RNGs/all.h"

/*
	currently the only non-uniform distribution supported is the Gaussian distribution
*/

namespace PractRand {
	namespace Internals {

		/*
			Gaussian/Normal distribution implementation

			Rolled my own.  
			Possibly a standard Ziggaraut method would work better, but this is simple and it works and it can be replaced later. 
			I suppose this does have the merit of using a fixed number of input bits.  

			At GAUSSIAN_CDF_TABLE_SIZE=256 this needs 2 KB of memory.  I'm not sure if that's too much or not enough.  
			It seems like a waste to spend 1 KB on something that many people will never use, but the linker may be smart enough to discard it.  
			If I went to a 16 KB table size then I could speed it up a tad without compromising quality too much, but this is fast enough for now.  

			I'm aiming for 25 or so bits of good resolution plus another 25 or so bits of noise here - 
			good enough to be practically impossible to distinguish from from true gaussian, and significantly better than single precision floats, but not quite maxing 
			out what double precision numbers can manage.  More than that and I would have have a hard time verifying things without large number types 
			and even more exotic math.  
		*/
		enum {GAUSSIAN_CDF_TABLE_SIZE_L2=8, GAUSSIAN_CDF_TABLE_SIZE = 1 << GAUSSIAN_CDF_TABLE_SIZE_L2};
		static const float primary_gaussian_cdf_table[GAUSSIAN_CDF_TABLE_SIZE] = {
			-1.6163263225, -1.3196592415, -1.1954910794, -1.1114453849, -1.0462748154, -0.9923525680, -0.9459735161, -0.9050338462, -0.8682159066, -0.8346362318, -0.8036718505, -0.7748660434, -0.7478734240, -0.7224260936, -0.6983118212, -0.6753594409, 
			-0.6534287703, -0.6324034597, -0.6121858041, -0.5926929007, -0.5738537537, -0.5556070563, -0.5378994703, -0.5206842715, -0.5039202731, -0.4875709602, -0.4716037875, -0.4559896044, -0.4407021805, -0.4257178113, -0.4110149884, -0.3965741208, 
			-0.3823772996, -0.3684080969, -0.3546513939, -0.3410932322, -0.3277206858, -0.3145217497, -0.3014852421, -0.2886007196, -0.2758584016, -0.2632491046, -0.2507641824, -0.2383954749, -0.2261352601, -0.2139762133, -0.2019113685, -0.1899340848, 
			-0.1780380153, -0.1662170791, -0.1544654357, -0.1427774615, -0.1311477279, -0.1195709818, -0.1080421265, -0.0965562049, -0.0851083829, -0.0736939343, -0.0623082262, -0.0509467056, -0.0396048858, -0.0282783338, -0.0169626579, -0.0056534962, 
			+0.0056534962, +0.0169626579, +0.0282783338, +0.0396048858, +0.0509467056, +0.0623082262, +0.0736939343, +0.0851083829, +0.0965562049, +0.1080421265, +0.1195709818, +0.1311477279, +0.1427774615, +0.1544654357, +0.1662170791, +0.1780380153, 
			+0.1899340848, +0.2019113685, +0.2139762133, +0.2261352601, +0.2383954749, +0.2507641824, +0.2632491046, +0.2758584016, +0.2886007196, +0.3014852421, +0.3145217497, +0.3277206858, +0.3410932322, +0.3546513939, +0.3684080969, +0.3823772996, 
			+0.3965741208, +0.4110149884, +0.4257178113, +0.4407021805, +0.4559896044, +0.4716037875, +0.4875709602, +0.5039202731, +0.5206842715, +0.5378994703, +0.5556070563, +0.5738537537, +0.5926929007, +0.6121858041, +0.6324034597, +0.6534287703, 
			+0.6753594409, +0.6983118212, +0.7224260936, +0.7478734240, +0.7748660434, +0.8036718505, +0.8346362318, +0.8682159066, +0.9050338462, +0.9459735161, +0.9923525680, +1.0462748154, +1.1114453849, +1.1954910794, +1.3196592415, +1.6163263225
		};
		static const float secondary_gaussian_cdf_table[GAUSSIAN_CDF_TABLE_SIZE] = {
			+1.0270065287e-010, +3.5445721984e-011, +2.2374584217e-011, +1.6762242066e-011, +1.3585169202e-011, +1.1524333934e-011, +1.0072595101e-011, +8.9914242845e-012, +8.1532648440e-012, +7.4834897566e-012, +6.9354238958e-012, +6.4783253196e-012, +6.0910925944e-012, +5.7587440258e-012, +5.4703391261e-012, +5.2176958089e-012,
			+4.9945685304e-012, +4.7961046655e-012, +4.6184748659e-012, +4.4586155873e-012, +4.3140458994e-012, +4.1827346681e-012, +4.0630026198e-012, +3.9534490258e-012, +3.8528960586e-012, +3.7603460307e-012, +3.6749481583e-012, +3.5959724575e-012, +3.5227890436e-012, +3.4548515714e-012, +3.3916838777e-012, +3.3328691249e-012,
			+3.2780409160e-012, +3.2268759731e-012, +3.1790880693e-012, +3.1344229687e-012, +3.0926541857e-012, +3.0535794117e-012, +3.0170174920e-012, +2.9828058540e-012, +2.9507983138e-012, +2.9208631954e-012, +2.8928817138e-012, +2.8667465801e-012, +2.8423607943e-012, +2.8196365979e-012, +2.7984945628e-012, +2.7788627969e-012,
			+2.7606762508e-012, +2.7438761111e-012, +2.7284092696e-012, +2.7142278580e-012, +2.7012888411e-012, +2.6895536597e-012, +2.6789879199e-012, +2.6695611208e-012, +2.6612464195e-012, +2.6540204272e-012, +2.6478630353e-012, +2.6427572677e-012, +2.6386891581e-012, +2.6356476497e-012, +2.6336245174e-012, +2.6326143094e-012,
			+2.6326143094e-012, +2.6336245174e-012, +2.6356476497e-012, +2.6386891581e-012, +2.6427572677e-012, +2.6478630353e-012, +2.6540204272e-012, +2.6612464195e-012, +2.6695611208e-012, +2.6789879199e-012, +2.6895536597e-012, +2.7012888411e-012, +2.7142278580e-012, +2.7284092696e-012, +2.7438761111e-012, +2.7606762508e-012,
			+2.7788627969e-012, +2.7984945628e-012, +2.8196365979e-012, +2.8423607943e-012, +2.8667465801e-012, +2.8928817138e-012, +2.9208631954e-012, +2.9507983138e-012, +2.9828058540e-012, +3.0170174920e-012, +3.0535794117e-012, +3.0926541857e-012, +3.1344229687e-012, +3.1790880693e-012, +3.2268759731e-012, +3.2780409160e-012,
			+3.3328691249e-012, +3.3916838777e-012, +3.4548515714e-012, +3.5227890436e-012, +3.5959724575e-012, +3.6749481583e-012, +3.7603460307e-012, +3.8528960586e-012, +3.9534490258e-012, +4.0630026198e-012, +4.1827346681e-012, +4.3140458994e-012, +4.4586155873e-012, +4.6184748659e-012, +4.7961046655e-012, +4.9945685304e-012,
			+5.2176958089e-012, +5.4703391261e-012, +5.7587440258e-012, +6.0910925944e-012, +6.4783253196e-012, +6.9354238958e-012, +7.4834897566e-012, +8.1532648440e-012, +8.9914242845e-012, +1.0072595101e-011, +1.1524333934e-011, +1.3585169202e-011, +1.6762242066e-011, +2.2374584217e-011, +3.5445721984e-011, +1.0270065287e-010
		};
		double generate_gaussian_fast(Uint64 raw64) {//fast CDF-based hybrid method
			Sint32 si = Sint32(raw64 >> 32);
			Uint32 indeces = Uint32(raw64);
			long index = (indeces >> (GAUSSIAN_CDF_TABLE_SIZE_L2*0)) & (GAUSSIAN_CDF_TABLE_SIZE-1);
			double rv = primary_gaussian_cdf_table[index];
			rv += si * secondary_gaussian_cdf_table[index];
			index = (indeces >> (GAUSSIAN_CDF_TABLE_SIZE_L2*1)) & (GAUSSIAN_CDF_TABLE_SIZE-1);
			rv += primary_gaussian_cdf_table[index];
			index = (indeces >> (GAUSSIAN_CDF_TABLE_SIZE_L2*2)) & (GAUSSIAN_CDF_TABLE_SIZE-1);
			rv += primary_gaussian_cdf_table[index];
			return rv;
		}
	}
}
