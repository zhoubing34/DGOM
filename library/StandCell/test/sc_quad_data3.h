#define NP 16
#define NFP 16

double quad_r[NP] = {-1.0000000000000000, -0.4472135954999578, 0.4472135954999578, 1.0000000000000000, -1.0000000000000000, -0.4472135954999578, 0.4472135954999578, 1.0000000000000000, -1.0000000000000000, -0.4472135954999578, 0.4472135954999578, 1.0000000000000000, -1.0000000000000000, -0.4472135954999578, 0.4472135954999578, 1.0000000000000000};
double quad_s[NP] = {-1.0000000000000000, -1.0000000000000000, -1.0000000000000000, -1.0000000000000000, -0.4472135954999578, -0.4472135954999578, -0.4472135954999578, -0.4472135954999578, 0.4472135954999578, 0.4472135954999578, 0.4472135954999578, 0.4472135954999578, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000};

double quad_V[NP][NP] = {{0.5000000000000001, -0.8660254037844386, 1.1180339887498949, -1.3228756555322954, -0.8660254037844386, 1.4999999999999998, -1.9364916731037083, 2.2912878474779199, 1.1180339887498949, -1.9364916731037083, 2.5000000000000004, -2.9580398915498081, -1.3228756555322954, 2.2912878474779199, -2.9580398915498081, 3.5000000000000000, },
                         {0.5000000000000001, -0.8660254037844386, 1.1180339887498949, -1.3228756555322954, -0.3872983346207416, 0.6708203932499367, -0.8660254037844385, 1.0246950765959595, -0.2236067977499792, 0.3872983346207420, -0.5000000000000004, 0.5916079783099621, 0.5916079783099617, -1.0246950765959597, 1.3228756555322954, -1.5652475842498528, },
                         {0.5000000000000001, -0.8660254037844386, 1.1180339887498949, -1.3228756555322954, 0.3872983346207416, -0.6708203932499367, 0.8660254037844385, -1.0246950765959595, -0.2236067977499792, 0.3872983346207420, -0.5000000000000004, 0.5916079783099621, -0.5916079783099617, 1.0246950765959597, -1.3228756555322954, 1.5652475842498528, },
                         {0.5000000000000001, -0.8660254037844386, 1.1180339887498949, -1.3228756555322954, 0.8660254037844386, -1.4999999999999998, 1.9364916731037083, -2.2912878474779199, 1.1180339887498949, -1.9364916731037083, 2.5000000000000004, -2.9580398915498081, 1.3228756555322954, -2.2912878474779199, 2.9580398915498081, -3.5000000000000000, },
                         {0.5000000000000001, -0.3872983346207416, -0.2236067977499792, 0.5916079783099617, -0.8660254037844386, 0.6708203932499367, 0.3872983346207420, -1.0246950765959597, 1.1180339887498949, -0.8660254037844385, -0.5000000000000004, 1.3228756555322954, -1.3228756555322954, 1.0246950765959595, 0.5916079783099621, -1.5652475842498528, },
                         {0.5000000000000001, -0.3872983346207416, -0.2236067977499792, 0.5916079783099617, -0.3872983346207416, 0.2999999999999998, 0.1732050807568878, -0.4582575694955839, -0.2236067977499792, 0.1732050807568878, 0.1000000000000002, -0.2645751311064593, 0.5916079783099617, -0.4582575694955839, -0.2645751311064593, 0.7000000000000001, },
                         {0.5000000000000001, -0.3872983346207416, -0.2236067977499792, 0.5916079783099617, 0.3872983346207416, -0.2999999999999998, -0.1732050807568878, 0.4582575694955839, -0.2236067977499792, 0.1732050807568878, 0.1000000000000002, -0.2645751311064593, -0.5916079783099617, 0.4582575694955839, 0.2645751311064593, -0.7000000000000001, },
                         {0.5000000000000001, -0.3872983346207416, -0.2236067977499792, 0.5916079783099617, 0.8660254037844386, -0.6708203932499367, -0.3872983346207420, 1.0246950765959597, 1.1180339887498949, -0.8660254037844385, -0.5000000000000004, 1.3228756555322954, 1.3228756555322954, -1.0246950765959595, -0.5916079783099621, 1.5652475842498528, },
                         {0.5000000000000001, 0.3872983346207416, -0.2236067977499792, -0.5916079783099617, -0.8660254037844386, -0.6708203932499367, 0.3872983346207420, 1.0246950765959597, 1.1180339887498949, 0.8660254037844385, -0.5000000000000004, -1.3228756555322954, -1.3228756555322954, -1.0246950765959595, 0.5916079783099621, 1.5652475842498528, },
                         {0.5000000000000001, 0.3872983346207416, -0.2236067977499792, -0.5916079783099617, -0.3872983346207416, -0.2999999999999998, 0.1732050807568878, 0.4582575694955839, -0.2236067977499792, -0.1732050807568878, 0.1000000000000002, 0.2645751311064593, 0.5916079783099617, 0.4582575694955839, -0.2645751311064593, -0.7000000000000001, },
                         {0.5000000000000001, 0.3872983346207416, -0.2236067977499792, -0.5916079783099617, 0.3872983346207416, 0.2999999999999998, -0.1732050807568878, -0.4582575694955839, -0.2236067977499792, -0.1732050807568878, 0.1000000000000002, 0.2645751311064593, -0.5916079783099617, -0.4582575694955839, 0.2645751311064593, 0.7000000000000001, },
                         {0.5000000000000001, 0.3872983346207416, -0.2236067977499792, -0.5916079783099617, 0.8660254037844386, 0.6708203932499367, -0.3872983346207420, -1.0246950765959597, 1.1180339887498949, 0.8660254037844385, -0.5000000000000004, -1.3228756555322954, 1.3228756555322954, 1.0246950765959595, -0.5916079783099621, -1.5652475842498528, },
                         {0.5000000000000001, 0.8660254037844386, 1.1180339887498949, 1.3228756555322954, -0.8660254037844386, -1.4999999999999998, -1.9364916731037083, -2.2912878474779199, 1.1180339887498949, 1.9364916731037083, 2.5000000000000004, 2.9580398915498081, -1.3228756555322954, -2.2912878474779199, -2.9580398915498081, -3.5000000000000000, },
                         {0.5000000000000001, 0.8660254037844386, 1.1180339887498949, 1.3228756555322954, -0.3872983346207416, -0.6708203932499367, -0.8660254037844385, -1.0246950765959595, -0.2236067977499792, -0.3872983346207420, -0.5000000000000004, -0.5916079783099621, 0.5916079783099617, 1.0246950765959597, 1.3228756555322954, 1.5652475842498528, },
                         {0.5000000000000001, 0.8660254037844386, 1.1180339887498949, 1.3228756555322954, 0.3872983346207416, 0.6708203932499367, 0.8660254037844385, 1.0246950765959595, -0.2236067977499792, -0.3872983346207420, -0.5000000000000004, -0.5916079783099621, -0.5916079783099617, -1.0246950765959597, -1.3228756555322954, -1.5652475842498528, },
                         {0.5000000000000001, 0.8660254037844386, 1.1180339887498949, 1.3228756555322954, 0.8660254037844386, 1.4999999999999998, 1.9364916731037083, 2.2912878474779199, 1.1180339887498949, 1.9364916731037083, 2.5000000000000004, 2.9580398915498081, 1.3228756555322954, 2.2912878474779199, 2.9580398915498081, 3.5000000000000000}};

double quad_M[NP][NP] = {{0.0204081632653061, 0.0076056733928565, -0.0076056733928564, 0.0034013605442177, 0.0076056733928565, 0.0028344671201815, -0.0028344671201814, 0.0012676122321427, -0.0076056733928564, -0.0028344671201814, 0.0028344671201814, -0.0012676122321427, 0.0034013605442177, 0.0012676122321428, -0.0012676122321427, 0.0005668934240363, },
                         {0.0076056733928565, 0.1020408163265306, 0.0170068027210884, -0.0076056733928564, 0.0028344671201814, 0.0380283669642823, 0.0063380611607137, -0.0028344671201814, -0.0028344671201814, -0.0380283669642821, -0.0063380611607137, 0.0028344671201814, 0.0012676122321427, 0.0170068027210884, 0.0028344671201814, -0.0012676122321427, },
                         {-0.0076056733928564, 0.0170068027210884, 0.1020408163265306, 0.0076056733928564, -0.0028344671201814, 0.0063380611607137, 0.0380283669642822, 0.0028344671201814, 0.0028344671201814, -0.0063380611607137, -0.0380283669642821, -0.0028344671201814, -0.0012676122321427, 0.0028344671201814, 0.0170068027210884, 0.0012676122321427, },
                         {0.0034013605442177, -0.0076056733928564, 0.0076056733928564, 0.0204081632653061, 0.0012676122321427, -0.0028344671201814, 0.0028344671201814, 0.0076056733928564, -0.0012676122321427, 0.0028344671201814, -0.0028344671201814, -0.0076056733928564, 0.0005668934240363, -0.0012676122321427, 0.0012676122321427, 0.0034013605442177, },
                         {0.0076056733928565, 0.0028344671201814, -0.0028344671201814, 0.0012676122321427, 0.1020408163265307, 0.0380283669642823, -0.0380283669642820, 0.0170068027210884, 0.0170068027210885, 0.0063380611607139, -0.0063380611607136, 0.0028344671201814, -0.0076056733928564, -0.0028344671201814, 0.0028344671201814, -0.0012676122321427, },
                         {0.0028344671201815, 0.0380283669642823, 0.0063380611607137, -0.0028344671201814, 0.0380283669642823, 0.5102040816326529, 0.0850340136054420, -0.0380283669642821, 0.0063380611607136, 0.0850340136054421, 0.0141723356009067, -0.0063380611607137, -0.0028344671201815, -0.0380283669642821, -0.0063380611607138, 0.0028344671201814, },
                         {-0.0028344671201814, 0.0063380611607137, 0.0380283669642822, 0.0028344671201814, -0.0380283669642820, 0.0850340136054420, 0.5102040816326530, 0.0380283669642822, -0.0063380611607138, 0.0141723356009069, 0.0850340136054420, 0.0063380611607136, 0.0028344671201813, -0.0063380611607138, -0.0380283669642822, -0.0028344671201814, },
                         {0.0012676122321427, -0.0028344671201814, 0.0028344671201814, 0.0076056733928564, 0.0170068027210884, -0.0380283669642821, 0.0380283669642822, 0.1020408163265306, 0.0028344671201814, -0.0063380611607137, 0.0063380611607137, 0.0170068027210884, -0.0012676122321427, 0.0028344671201814, -0.0028344671201814, -0.0076056733928564, },
                         {-0.0076056733928564, -0.0028344671201814, 0.0028344671201814, -0.0012676122321427, 0.0170068027210885, 0.0063380611607136, -0.0063380611607138, 0.0028344671201814, 0.1020408163265306, 0.0380283669642821, -0.0380283669642822, 0.0170068027210884, 0.0076056733928564, 0.0028344671201814, -0.0028344671201814, 0.0012676122321428, },
                         {-0.0028344671201814, -0.0380283669642821, -0.0063380611607137, 0.0028344671201814, 0.0063380611607139, 0.0850340136054421, 0.0141723356009069, -0.0063380611607137, 0.0380283669642821, 0.5102040816326532, 0.0850340136054421, -0.0380283669642821, 0.0028344671201813, 0.0380283669642823, 0.0063380611607136, -0.0028344671201813, },
                         {0.0028344671201814, -0.0063380611607137, -0.0380283669642821, -0.0028344671201814, -0.0063380611607136, 0.0141723356009067, 0.0850340136054420, 0.0063380611607137, -0.0380283669642822, 0.0850340136054421, 0.5102040816326529, 0.0380283669642822, -0.0028344671201815, 0.0063380611607137, 0.0380283669642822, 0.0028344671201815, },
                         {-0.0012676122321427, 0.0028344671201814, -0.0028344671201814, -0.0076056733928564, 0.0028344671201814, -0.0063380611607137, 0.0063380611607136, 0.0170068027210884, 0.0170068027210884, -0.0380283669642821, 0.0380283669642822, 0.1020408163265306, 0.0012676122321427, -0.0028344671201814, 0.0028344671201814, 0.0076056733928564, },
                         {0.0034013605442177, 0.0012676122321427, -0.0012676122321427, 0.0005668934240363, -0.0076056733928564, -0.0028344671201815, 0.0028344671201813, -0.0012676122321427, 0.0076056733928564, 0.0028344671201813, -0.0028344671201815, 0.0012676122321427, 0.0204081632653061, 0.0076056733928564, -0.0076056733928564, 0.0034013605442177, },
                         {0.0012676122321428, 0.0170068027210884, 0.0028344671201814, -0.0012676122321427, -0.0028344671201814, -0.0380283669642821, -0.0063380611607138, 0.0028344671201814, 0.0028344671201814, 0.0380283669642823, 0.0063380611607137, -0.0028344671201814, 0.0076056733928564, 0.1020408163265306, 0.0170068027210884, -0.0076056733928564, },
                         {-0.0012676122321427, 0.0028344671201814, 0.0170068027210884, 0.0012676122321427, 0.0028344671201814, -0.0063380611607138, -0.0380283669642822, -0.0028344671201814, -0.0028344671201814, 0.0063380611607136, 0.0380283669642822, 0.0028344671201814, -0.0076056733928564, 0.0170068027210884, 0.1020408163265306, 0.0076056733928564, },
                         {0.0005668934240363, -0.0012676122321427, 0.0012676122321427, 0.0034013605442177, -0.0012676122321427, 0.0028344671201814, -0.0028344671201814, -0.0076056733928564, 0.0012676122321428, -0.0028344671201813, 0.0028344671201815, 0.0076056733928564, 0.0034013605442177, -0.0076056733928564, 0.0076056733928564, 0.0204081632653061, }};

double quad_Dr[NP][NP] =  {{-2.9999999999999991, 4.0450849718747373, -1.5450849718747377, 0.4999999999999999, -0.0000000000000002, -0.0000000000000004, 0.0000000000000008, -0.0000000000000002, -0.0000000000000001, 0.0000000000000009, -0.0000000000000009, 0.0000000000000001, 0.0000000000000000, 0.0000000000000001, -0.0000000000000003, 0.0000000000000001, },
                           {-0.8090169943749471, -0.0000000000000009, 1.1180339887498951, -0.3090169943749471, 0.0000000000000000, 0.0000000000000003, 0.0000000000000001, -0.0000000000000003, 0.0000000000000000, -0.0000000000000003, 0.0000000000000003, -0.0000000000000000, -0.0000000000000000, 0.0000000000000002, -0.0000000000000003, 0.0000000000000001, },
                           {0.3090169943749475, -1.1180339887498951, 0.0000000000000009, 0.8090169943749470, 0.0000000000000001, -0.0000000000000004, 0.0000000000000000, 0.0000000000000002, -0.0000000000000000, -0.0000000000000001, 0.0000000000000001, 0.0000000000000000, -0.0000000000000002, 0.0000000000000001, 0.0000000000000001, 0.0000000000000000, },
                           {-0.4999999999999993, 1.5450849718747379, -4.0450849718747381, 3.0000000000000000, -0.0000000000000002, -0.0000000000000012, 0.0000000000000024, -0.0000000000000011, -0.0000000000000002, 0.0000000000000014, -0.0000000000000014, 0.0000000000000002, 0.0000000000000001, -0.0000000000000007, 0.0000000000000008, -0.0000000000000002, },
                           {0.0000000000000012, 0.0000000000000004, -0.0000000000000005, 0.0000000000000002, -2.9999999999999991, 4.0450849718747364, -1.5450849718747370, 0.4999999999999999, -0.0000000000000009, 0.0000000000000018, -0.0000000000000018, 0.0000000000000010, 0.0000000000000001, -0.0000000000000003, 0.0000000000000004, -0.0000000000000003, },
                           {0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000001, -0.8090169943749471, -0.0000000000000007, 1.1180339887498951, -0.3090169943749472, -0.0000000000000001, 0.0000000000000001, 0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000000, },
                           {0.0000000000000000, -0.0000000000000000, 0.0000000000000001, -0.0000000000000001, 0.3090169943749473, -1.1180339887498951, 0.0000000000000008, 0.8090169943749470, 0.0000000000000001, -0.0000000000000002, 0.0000000000000001, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, -0.0000000000000000, 0.0000000000000001, },
                           {0.0000000000000000, 0.0000000000000002, -0.0000000000000000, -0.0000000000000002, -0.4999999999999996, 1.5450849718747373, -4.0450849718747364, 2.9999999999999987, -0.0000000000000004, 0.0000000000000004, -0.0000000000000004, 0.0000000000000003, -0.0000000000000000, -0.0000000000000000, -0.0000000000000001, 0.0000000000000001, },
                           {0.0000000000000000, 0.0000000000000000, 0.0000000000000002, -0.0000000000000000, 0.0000000000000000, 0.0000000000000013, -0.0000000000000011, 0.0000000000000002, -3.0000000000000000, 4.0450849718747373, -1.5450849718747377, 0.5000000000000003, 0.0000000000000000, -0.0000000000000005, 0.0000000000000000, -0.0000000000000002, },
                           {0.0000000000000001, 0.0000000000000001, -0.0000000000000001, 0.0000000000000000, -0.0000000000000004, 0.0000000000000001, 0.0000000000000003, -0.0000000000000000, -0.8090169943749469, -0.0000000000000008, 1.1180339887498951, -0.3090169943749473, -0.0000000000000001, 0.0000000000000001, -0.0000000000000000, 0.0000000000000000, },
                           {0.0000000000000000, 0.0000000000000000, 0.0000000000000001, 0.0000000000000000, 0.0000000000000001, -0.0000000000000002, -0.0000000000000001, 0.0000000000000001, 0.3090169943749473, -1.1180339887498953, 0.0000000000000011, 0.8090169943749469, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.0000000000000000, },
                           {0.0000000000000002, -0.0000000000000004, 0.0000000000000004, 0.0000000000000002, 0.0000000000000006, 0.0000000000000005, -0.0000000000000011, -0.0000000000000001, -0.5000000000000004, 1.5450849718747373, -4.0450849718747364, 2.9999999999999996, -0.0000000000000007, 0.0000000000000006, -0.0000000000000002, -0.0000000000000001, },
                           {0.0000000000000004, -0.0000000000000002, 0.0000000000000002, -0.0000000000000002, -0.0000000000000003, 0.0000000000000006, 0.0000000000000009, -0.0000000000000012, -0.0000000000000002, 0.0000000000000006, -0.0000000000000007, 0.0000000000000003, -3.0000000000000000, 4.0450849718747364, -1.5450849718747373, 0.5000000000000003, },
                           {-0.0000000000000002, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -0.0000000000000002, 0.0000000000000001, -0.0000000000000000, 0.0000000000000002, 0.0000000000000000, 0.0000000000000001, 0.0000000000000000, -0.0000000000000001, -0.8090169943749471, -0.0000000000000007, 1.1180339887498951, -0.3090169943749473, },
                           {0.0000000000000001, 0.0000000000000000, -0.0000000000000001, 0.0000000000000000, -0.0000000000000002, -0.0000000000000001, 0.0000000000000002, 0.0000000000000001, 0.0000000000000002, -0.0000000000000003, 0.0000000000000002, -0.0000000000000001, 0.3090169943749472, -1.1180339887498951, 0.0000000000000007, 0.8090169943749471, },
                           {-0.0000000000000001, 0.0000000000000000, 0.0000000000000001, -0.0000000000000003, -0.0000000000000001, 0.0000000000000003, -0.0000000000000007, 0.0000000000000004, -0.0000000000000003, 0.0000000000000003, -0.0000000000000001, 0.0000000000000002, -0.5000000000000008, 1.5450849718747377, -4.0450849718747373, 3.0000000000000004, }};

double quad_Ds[NP][NP] = {{-2.9999999999999996, -0.0000000000000001, 0.0000000000000003, -0.0000000000000001, 4.0450849718747364, 0.0000000000000003, 0.0000000000000001, 0.0000000000000000, -1.5450849718747373, -0.0000000000000001, -0.0000000000000004, -0.0000000000000001, 0.5000000000000004, -0.0000000000000001, -0.0000000000000001, 0.0000000000000002, },
                          {0.0000000000000003, -3.0000000000000000, 0.0000000000000009, -0.0000000000000007, 0.0000000000000004, 4.0450849718747373, -0.0000000000000010, 0.0000000000000008, -0.0000000000000001, -1.5450849718747377, 0.0000000000000000, -0.0000000000000001, -0.0000000000000004, 0.5000000000000001, 0.0000000000000004, -0.0000000000000002, },
                          {0.0000000000000002, -0.0000000000000004, -2.9999999999999987, -0.0000000000000008, -0.0000000000000001, 0.0000000000000002, 4.0450849718747364, 0.0000000000000012, -0.0000000000000000, 0.0000000000000005, -1.5450849718747381, -0.0000000000000002, -0.0000000000000002, -0.0000000000000002, 0.5000000000000008, -0.0000000000000002, },
                          {-0.0000000000000001, 0.0000000000000002, 0.0000000000000000, -2.9999999999999996, 0.0000000000000002, 0.0000000000000003, 0.0000000000000001, 4.0450849718747364, -0.0000000000000002, 0.0000000000000000, -0.0000000000000005, -1.5450849718747373, 0.0000000000000001, -0.0000000000000005, 0.0000000000000003, 0.5000000000000004, },
                          {-0.8090169943749472, -0.0000000000000001, -0.0000000000000000, 0.0000000000000000, -0.0000000000000007, 0.0000000000000001, 0.0000000000000000, -0.0000000000000000, 1.1180339887498951, 0.0000000000000001, 0.0000000000000000, -0.0000000000000000, -0.3090169943749472, -0.0000000000000001, -0.0000000000000000, 0.0000000000000000, },
                          {-0.0000000000000001, -0.8090169943749472, 0.0000000000000002, -0.0000000000000002, -0.0000000000000001, -0.0000000000000007, -0.0000000000000001, 0.0000000000000002, 0.0000000000000001, 1.1180339887498951, -0.0000000000000002, 0.0000000000000001, -0.0000000000000001, -0.3090169943749470, -0.0000000000000002, 0.0000000000000001, },
                          {-0.0000000000000002, 0.0000000000000002, -0.8090169943749471, -0.0000000000000001, 0.0000000000000001, -0.0000000000000004, -0.0000000000000004, 0.0000000000000000, 0.0000000000000001, 0.0000000000000003, 1.1180339887498947, 0.0000000000000001, -0.0000000000000001, 0.0000000000000001, -0.3090169943749473, -0.0000000000000001, },
                          {-0.0000000000000002, -0.0000000000000005, 0.0000000000000004, -0.8090169943749470, 0.0000000000000002, 0.0000000000000006, -0.0000000000000006, -0.0000000000000008, 0.0000000000000000, -0.0000000000000001, 0.0000000000000002, 1.1180339887498951, 0.0000000000000000, 0.0000000000000000, -0.0000000000000001, -0.3090169943749473, },
                          {0.3090169943749472, 0.0000000000000001, -0.0000000000000000, 0.0000000000000000, -1.1180339887498951, -0.0000000000000001, -0.0000000000000002, 0.0000000000000000, 0.0000000000000007, -0.0000000000000000, 0.0000000000000002, 0.0000000000000000, 0.8090169943749472, 0.0000000000000001, 0.0000000000000000, -0.0000000000000001, },
                          {0.0000000000000002, 0.3090169943749472, 0.0000000000000000, 0.0000000000000000, -0.0000000000000002, -1.1180339887498951, 0.0000000000000001, -0.0000000000000001, 0.0000000000000000, 0.0000000000000009, -0.0000000000000000, -0.0000000000000001, 0.0000000000000000, 0.8090169943749469, 0.0000000000000001, 0.0000000000000001, },
                          {-0.0000000000000000, -0.0000000000000001, 0.3090169943749474, 0.0000000000000001, -0.0000000000000001, 0.0000000000000000, -1.1180339887498951, -0.0000000000000001, 0.0000000000000001, -0.0000000000000000, 0.0000000000000009, -0.0000000000000001, 0.0000000000000002, 0.0000000000000000, 0.8090169943749470, 0.0000000000000000, },
                          {0.0000000000000002, -0.0000000000000001, 0.0000000000000001, 0.3090169943749472, -0.0000000000000002, -0.0000000000000000, -0.0000000000000002, -1.1180339887498949, 0.0000000000000000, 0.0000000000000002, -0.0000000000000001, 0.0000000000000007, 0.0000000000000001, -0.0000000000000001, 0.0000000000000001, 0.8090169943749470, },
                          {-0.4999999999999996, -0.0000000000000002, -0.0000000000000000, -0.0000000000000000, 1.5450849718747368, 0.0000000000000006, 0.0000000000000002, -0.0000000000000002, -4.0450849718747364, -0.0000000000000007, -0.0000000000000000, 0.0000000000000002, 2.9999999999999991, 0.0000000000000003, -0.0000000000000001, -0.0000000000000001, },
                          {-0.0000000000000003, -0.5000000000000004, 0.0000000000000000, -0.0000000000000001, 0.0000000000000002, 1.5450849718747379, -0.0000000000000009, 0.0000000000000007, 0.0000000000000000, -4.0450849718747373, 0.0000000000000000, -0.0000000000000000, -0.0000000000000002, 3.0000000000000000, 0.0000000000000004, -0.0000000000000003, },
                          {-0.0000000000000000, 0.0000000000000000, -0.5000000000000002, 0.0000000000000003, -0.0000000000000007, 0.0000000000000011, 1.5450849718747361, 0.0000000000000005, 0.0000000000000001, -0.0000000000000014, -4.0450849718747355, -0.0000000000000001, 0.0000000000000006, -0.0000000000000001, 3.0000000000000000, -0.0000000000000007, },
                          {-0.0000000000000004, 0.0000000000000010, -0.0000000000000012, -0.4999999999999998, 0.0000000000000004, -0.0000000000000007, 0.0000000000000014, 1.5450849718747368, 0.0000000000000000, 0.0000000000000001, -0.0000000000000008, -4.0450849718747364, -0.0000000000000001, -0.0000000000000004, 0.0000000000000006, 2.9999999999999996, }};

double quad_LIFT[NP][NFP] = {{7.9999999999999973, -0.0000000000000119, -0.0000000000000040, 0.0000000000000000, -1.9999999999999989, 0.0000000000000033, 0.0000000000000017, 0.0000000000000008, 0.0000000000000010, -0.0000000000000011, -0.0000000000000011, -1.9999999999999989, -0.0000000000000004, 0.0000000000000053, 0.0000000000000000, 7.9999999999999973, },
                             {-0.0000000000000014, 8.0000000000000036, -0.0000000000000008, -0.0000000000000002, 0.8944271909999155, -0.0000000000000013, -0.0000000000000000, 0.0000000000000001, 0.0000000000000002, 0.0000000000000009, -1.9999999999999991, -0.0000000000000002, -0.0000000000000002, -0.0000000000000011, 0.0000000000000001, -0.8944271909999163, },
                             {-0.0000000000000007, -0.0000000000000037, 8.0000000000000000, 0.0000000000000002, -0.8944271909999163, 0.0000000000000014, -0.0000000000000001, -0.0000000000000002, 0.0000000000000002, -2.0000000000000009, 0.0000000000000003, 0.0000000000000000, 0.0000000000000001, 0.0000000000000011, 0.0000000000000002, 0.8944271909999154, },
                             {-0.0000000000000004, 0.0000000000000062, 0.0000000000000009, 7.9999999999999973, 7.9999999999999973, -0.0000000000000129, -0.0000000000000049, 0.0000000000000000, -1.9999999999999991, 0.0000000000000022, 0.0000000000000016, 0.0000000000000009, 0.0000000000000009, -0.0000000000000009, 0.0000000000000003, -1.9999999999999989, },
                             {-0.8944271909999163, 0.0000000000000013, 0.0000000000000006, 0.0000000000000001, 0.0000000000000004, -2.0000000000000009, 0.0000000000000006, 0.0000000000000001, 0.0000000000000001, 0.0000000000000010, -0.0000000000000004, 0.8944271909999154, -0.0000000000000007, -0.0000000000000038, 8.0000000000000000, 0.0000000000000001, },
                             {0.0000000000000001, -0.8944271909999173, 0.0000000000000002, 0.0000000000000001, -0.0000000000000001, 0.8944271909999165, 0.0000000000000001, 0.0000000000000001, 0.0000000000000000, -0.0000000000000002, 0.8944271909999157, 0.0000000000000000, 0.0000000000000001, 0.0000000000000005, -0.8944271909999166, -0.0000000000000000, },
                             {0.0000000000000001, 0.0000000000000004, -0.8944271909999167, -0.0000000000000000, 0.0000000000000001, -0.8944271909999173, 0.0000000000000002, 0.0000000000000000, -0.0000000000000001, 0.8944271909999165, 0.0000000000000000, 0.0000000000000001, -0.0000000000000000, -0.0000000000000002, 0.8944271909999160, 0.0000000000000001, },
                             {0.0000000000000002, -0.0000000000000006, 0.0000000000000000, -0.8944271909999162, -0.0000000000000014, 8.0000000000000036, -0.0000000000000013, -0.0000000000000002, 0.8944271909999154, -0.0000000000000018, -0.0000000000000001, 0.0000000000000001, 0.0000000000000002, 0.0000000000000012, -1.9999999999999991, 0.0000000000000001, },
                             {0.8944271909999155, -0.0000000000000018, -0.0000000000000001, 0.0000000000000001, 0.0000000000000003, 0.0000000000000012, -1.9999999999999989, 0.0000000000000001, 0.0000000000000002, -0.0000000000000006, -0.0000000000000000, -0.8944271909999162, -0.0000000000000016, 8.0000000000000053, -0.0000000000000013, -0.0000000000000001, },
                             {-0.0000000000000001, 0.8944271909999165, 0.0000000000000000, 0.0000000000000001, -0.0000000000000000, -0.0000000000000002, 0.8944271909999160, 0.0000000000000001, 0.0000000000000001, 0.0000000000000004, -0.8944271909999167, -0.0000000000000000, 0.0000000000000001, -0.8944271909999173, 0.0000000000000002, 0.0000000000000000, },
                             {-0.0000000000000000, -0.0000000000000003, 0.8944271909999157, 0.0000000000000000, 0.0000000000000001, 0.0000000000000005, -0.8944271909999167, -0.0000000000000000, 0.0000000000000001, -0.8944271909999174, 0.0000000000000002, 0.0000000000000001, -0.0000000000000001, 0.8944271909999165, 0.0000000000000000, 0.0000000000000001, },
                             {0.0000000000000000, 0.0000000000000011, -0.0000000000000004, 0.8944271909999154, -0.0000000000000007, -0.0000000000000037, 8.0000000000000000, 0.0000000000000002, -0.8944271909999163, 0.0000000000000013, 0.0000000000000006, 0.0000000000000001, 0.0000000000000004, -2.0000000000000009, 0.0000000000000005, 0.0000000000000001, },
                             {-1.9999999999999993, 0.0000000000000022, 0.0000000000000015, 0.0000000000000009, 0.0000000000000009, -0.0000000000000009, 0.0000000000000002, -1.9999999999999989, -0.0000000000000004, 0.0000000000000062, 0.0000000000000010, 7.9999999999999973, 7.9999999999999982, -0.0000000000000129, -0.0000000000000053, -0.0000000000000002, },
                             {0.0000000000000002, -2.0000000000000009, 0.0000000000000004, 0.0000000000000001, 0.0000000000000000, 0.0000000000000011, 0.0000000000000002, 0.8944271909999155, -0.0000000000000007, -0.0000000000000038, 8.0000000000000000, 0.0000000000000001, -0.8944271909999163, 0.0000000000000014, -0.0000000000000001, -0.0000000000000002, },
                             {0.0000000000000002, 0.0000000000000009, -1.9999999999999991, -0.0000000000000002, -0.0000000000000002, -0.0000000000000011, 0.0000000000000001, -0.8944271909999162, -0.0000000000000016, 8.0000000000000053, -0.0000000000000009, -0.0000000000000001, 0.8944271909999155, -0.0000000000000013, -0.0000000000000001, 0.0000000000000001, },
                             {0.0000000000000011, -0.0000000000000011, -0.0000000000000011, -1.9999999999999991, -0.0000000000000004, 0.0000000000000053, 0.0000000000000000, 7.9999999999999982, 7.9999999999999982, -0.0000000000000120, -0.0000000000000036, -0.0000000000000002, -1.9999999999999991, 0.0000000000000033, 0.0000000000000018, 0.0000000000000008, }};
