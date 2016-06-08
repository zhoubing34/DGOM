#ifndef QUAD03
#define QUAD03

double p_r[16] = {-1.00000000000000e+00,-4.47213595499958e-01,4.47213595499958e-01,1.00000000000000e+00,-1.00000000000000e+00,-4.47213595499958e-01,4.47213595499958e-01,1.00000000000000e+00,-1.00000000000000e+00,-4.47213595499958e-01,4.47213595499958e-01,1.00000000000000e+00,-1.00000000000000e+00,-4.47213595499958e-01,4.47213595499958e-01,1.00000000000000e+00};
double p_s[16] = {-1.00000000000000e+00,-1.00000000000000e+00,-1.00000000000000e+00,-1.00000000000000e+00,-4.47213595499958e-01,-4.47213595499958e-01,-4.47213595499958e-01,-4.47213595499958e-01,4.47213595499958e-01,4.47213595499958e-01,4.47213595499958e-01,4.47213595499958e-01,1.00000000000000e+00,1.00000000000000e+00,1.00000000000000e+00,1.00000000000000e+00};

double p_w[4] = {1.66666666666667e-01,8.33333333333333e-01,8.33333333333333e-01,1.66666666666667e-01};

double p_we[16] = {2.77777777777781e-02,1.38888888888889e-01,1.38888888888889e-01,2.77777777777777e-02,1.38888888888890e-01,6.94444444444444e-01,6.94444444444444e-01,1.38888888888889e-01,1.38888888888889e-01,6.94444444444445e-01,6.94444444444444e-01,1.38888888888889e-01,2.77777777777774e-02,1.38888888888889e-01,1.38888888888889e-01,2.77777777777781e-02};
double p_Dr[16][16] = {{-3.00000000000000e+00,4.04508497187474e+00,-1.54508497187474e+00,5.00000000000000e-01,-1.55230317610889e-16,-3.72125619086060e-16,7.77156117237610e-16,-2.49800180540660e-16,-1.39662868479567e-16,9.31625432986270e-16,-9.31625432986270e-16,1.39662868479567e-16,4.56956992688404e-17,1.11022302462516e-16,-2.95495879809501e-16,1.38777878078145e-16},
{-8.09016994374947e-01,-9.19369056662435e-16,1.11803398874990e+00,-3.09016994374947e-01,7.16572514891872e-18,2.84267818815185e-16,5.55111512312578e-17,-3.46944695195361e-16,5.58415400179150e-19,-2.72115337474691e-16,2.72115337474691e-16,-5.58415400179150e-19,-4.10337677138216e-17,1.66533453693773e-16,-2.50399776250282e-16,1.24900090270330e-16},{3.09016994374948e-01,-1.11803398874990e+00,9.13032195979069e-16,8.09016994374947e-01,9.27450139013422e-17,-3.64050127167868e-16,4.92605083414944e-17,2.22044604925031e-16,-5.43003948669704e-18,-5.55747332572895e-17,5.55747332572895e-17,5.43003948669704e-18,-1.68213697995899e-16,9.23555870473038e-17,7.58581109485954e-17,0.00000000000000e+00},{-4.99999999999999e-01,1.54508497187474e+00,-4.04508497187474e+00,3.00000000000000e+00,-1.65848242998315e-16,-1.16641938655187e-15,2.44249065417534e-15,-1.11022302462516e-15,-1.66623866373069e-16,1.44616929350198e-15,-1.44616929350198e-15,1.66623866373069e-16,5.57642045683040e-17,-6.66133814775094e-16,8.32414215131821e-16,-2.22044604925031e-16},{1.22124532708767e-15,3.97205464519563e-16,-5.07535837438523e-16,2.02517878193087e-16,-3.00000000000000e+00,4.04508497187474e+00,-1.54508497187474e+00,5.00000000000000e-01,-9.38747565313801e-16,1.77635683940025e-15,-1.80940053252066e-15,9.71791258434211e-16,1.11022302462516e-16,-2.86875091600604e-16,3.97205464519563e-16,-3.22817254573205e-16},{2.77555756156289e-17,-2.76340858942944e-17,-2.01423061606973e-17,5.45741018951358e-17,-8.09016994374947e-01,-7.39673998990353e-16,1.11803398874990e+00,-3.09016994374947e-01,-6.59534976526788e-17,8.74204513243470e-17,1.94108090302955e-17,-4.08777627019636e-17,-2.77555756156289e-17,2.30708344886404e-17,2.47055575663513e-17,-4.18511328885158e-17},{4.16333634234434e-17,-4.96506830649454e-17,9.74270751199372e-17,-7.91945321096747e-17,3.09016994374947e-01,-1.11803398874990e+00,7.72965075129737e-16,8.09016994374947e-01,1.26015520512144e-16,-2.22044604925031e-16,1.15213344570389e-16,-1.91842601575019e-17,-4.16333634234434e-17,1.87429100995367e-18,-4.96506830649454e-17,7.92689508858240e-17},{2.77555756156289e-17,1.54812517361149e-16,-4.44821444421893e-17,-1.88046668226118e-16,-5.00000000000000e-01,1.54508497187474e+00,-4.04508497187474e+00,3.00000000000000e+00,-3.73465454825148e-16,4.44089209850063e-16,-4.11045516729653e-16,3.40421761704738e-16,-2.77555756156289e-17,-4.30982853550767e-17,-6.72320875638826e-17,8.32962949611983e-17},{0.00000000000000e+00,0.00000000000000e+00,2.22044604925031e-16,-2.77555756156289e-17,0.00000000000000e+00,1.26882631385732e-15,-1.08401421422546e-15,2.04890715817255e-16,-3.00000000000000e+00,4.04508497187474e+00,-1.54508497187474e+00,5.00000000000000e-01,0.00000000000000e+00,-4.70298020249759e-16,0.00000000000000e+00,-1.77135140201626e-16},{9.93013661298909e-17,6.98562208778555e-17,-1.11022302462516e-16,2.77555756156289e-17,-4.44089209850063e-16,6.19447562791828e-17,2.64451350956441e-16,-3.83572621234814e-17,-8.09016994374947e-01,-8.19104519051478e-16,1.11803398874990e+00,-3.09016994374947e-01,-9.93013661298909e-17,6.86155564311058e-17,-3.04867428227798e-17,1.06016865078525e-17},{3.06858096987851e-17,0.00000000000000e+00,5.24176207993920e-17,0.00000000000000e+00,1.11022302462516e-16,-1.58603289232165e-16,-5.24176207993920e-17,1.11022302462516e-16,3.09016994374947e-01,-1.11803398874990e+00,1.11022302462516e-15,8.09016994374947e-01,-3.06858096987851e-17,-5.54667823983524e-32,0.00000000000000e+00,0.00000000000000e+00},{1.53429048493925e-16,-4.44089209850063e-16,4.44089209850063e-16,2.22044604925031e-16,5.55111512312578e-16,4.75809867696496e-16,-1.05780540382576e-15,-8.48134920628197e-17,-5.00000000000000e-01,1.54508497187474e+00,-4.04508497187474e+00,3.00000000000000e+00,-7.08540560806504e-16,6.13716193975702e-16,-2.22044604925031e-16,-1.37231112862212e-16},{4.44089209850063e-16,-2.22044604925031e-16,1.59490624061064e-16,-1.64469113640144e-16,-2.81944057743921e-16,5.59499813900210e-16,8.88178419700125e-16,-1.16573417585641e-15,-2.03503416882307e-16,5.59499813900210e-16,-6.99374767375262e-16,3.43378370357359e-16,-3.00000000000000e+00,4.04508497187474e+00,-1.54508497187474e+00,5.00000000000000e-01},{-2.22044604925031e-16,0.00000000000000e+00,3.12769904319838e-17,1.86698609813205e-17,-2.48758343951361e-16,6.99374767375263e-17,-4.32237377111964e-17,2.22044604925031e-16,2.98950648605870e-17,6.99374767375262e-17,1.23259516440783e-32,-9.98325415981133e-17,-8.09016994374947e-01,-6.54187067495881e-16,1.11803398874990e+00,-3.09016994374947e-01},{5.55111512312578e-17,2.62232450059618e-17,-5.75002354379456e-17,2.90177753352860e-17,-2.01410971262090e-16,-1.31655936125457e-16,2.22044604925031e-16,1.11022302462516e-16,1.73082606375664e-16,-3.17091439671373e-16,2.47153962933846e-16,-1.03145129638138e-16,3.09016994374947e-01,-1.11803398874990e+00,7.34724515192532e-16,8.09016994374947e-01},{-1.11022302462516e-16,0.00000000000000e+00,6.25539808639676e-17,-2.59995791323891e-16,-5.77053020250737e-17,2.79749906950105e-16,-6.66133814775094e-16,4.44089209850063e-16,-3.09015743237962e-16,2.79749906950105e-16,-1.39874953475052e-16,1.69140789762910e-16,-5.00000000000001e-01,1.54508497187474e+00,-4.04508497187474e+00,3.00000000000000e+00}};

double p_Ds[16][16] = {{-3.00000000000000e+00,-9.58722832435283e-17,2.89131848755571e-16,-1.37422873064925e-16,4.04508497187474e+00,3.14593151678196e-16,1.17548374108805e-16,2.35868743233960e-17,-1.54508497187474e+00,-8.02977091063623e-17,-3.51843816680639e-16,-1.28366981605200e-16,5.00000000000000e-01,-1.38423159328306e-16,-5.48364061837368e-17,2.42202980346729e-16},
{2.51487495847620e-16,-3.00000000000000e+00,8.88178419700125e-16,-6.83855769365058e-16,4.33927351873010e-16,4.04508497187474e+00,-9.99200722162641e-16,7.87317975214662e-16,-1.29999488968315e-16,-1.54508497187474e+00,0.00000000000000e+00,-9.20451159567162e-17,-3.62509798310136e-16,5.00000000000000e-01,4.44089209850063e-16,-2.04322650335067e-16},{2.42126423462405e-16,-4.44089209850063e-16,-3.00000000000000e+00,-7.62075126777109e-16,-5.72957828880136e-17,1.91084094455466e-16,4.04508497187474e+00,1.16751880751317e-15,-1.29344807471578e-17,5.08290672919796e-16,-1.54508497187474e+00,-2.09110124177874e-16,-1.71896159827234e-16,-2.22044604925031e-16,5.00000000000001e-01,-1.96333556558187e-16},{-9.93013661298909e-17,1.77587390981351e-16,1.56721745306915e-17,-3.00000000000000e+00,2.22044604925031e-16,3.18286653154623e-16,1.13854872632379e-16,4.04508497187474e+00,-2.22044604925031e-16,2.50944721896054e-17,-4.57235997976607e-16,-1.54508497187474e+00,9.93013661298909e-17,-5.20968516325580e-16,3.27708950813537e-16,5.00000000000000e-01},{-8.09016994374947e-01,-5.83432055719806e-17,-3.33585135163924e-18,1.18625563597354e-17,-6.66133814775094e-16,5.83432055719805e-17,3.33585135163920e-18,-1.18625563597354e-17,1.11803398874990e+00,5.83432055719806e-17,3.33585135163927e-18,-1.18625563597354e-17,-3.09016994374947e-01,-5.83432055719805e-17,-3.33585135163922e-18,1.18625563597354e-17},{-7.62828979127028e-17,-8.09016994374947e-01,2.33765541257656e-16,-1.51622175178641e-16,-9.58711608294410e-17,-7.13714801544744e-16,-5.55111512312578e-17,1.51382312060699e-16,1.46295226502353e-16,1.11803398874990e+00,-2.22044604925031e-16,7.57493784226785e-17,-1.25331716604590e-16,-3.09016994374947e-01,-2.22044604925031e-16,7.56810335396446e-17},{-1.82844410350754e-16,1.66533453693773e-16,-8.09016994374947e-01,-5.09211309069017e-17,1.43872741808810e-16,-3.52175316833094e-16,-3.88578058618805e-16,2.26607118849637e-17,8.76542997993017e-17,2.82237840095567e-16,1.11803398874989e+00,1.34390305125730e-16,-9.24728461559902e-17,5.55111512312578e-17,-3.09016994374947e-01,-6.23396712051585e-17},{-1.84114858192711e-16,-5.00349656305187e-16,4.38670599381567e-16,-8.09016994374947e-01,1.66533453693773e-16,6.38283727591060e-16,-5.76604670667440e-16,-8.32667268468867e-16,0.00000000000000e+00,-1.40469811557416e-16,2.02148868481036e-16,1.11803398874990e+00,1.75814044989372e-17,2.53574027154283e-18,-6.42147971951626e-17,-3.09016994374947e-01},{3.09016994374947e-01,5.24742821496320e-17,-1.09851239723193e-17,2.65703547124383e-17,-1.11803398874990e+00,-1.03714089284542e-16,-1.51700941858456e-16,1.81087040477629e-17,7.00128461923931e-16,-2.79024443766036e-17,1.59959361672362e-16,1.22933624438732e-17,8.09016994374947e-01,7.91422515115141e-17,2.72670415841281e-18,-5.69724212040743e-17},{2.43711640070007e-16,3.09016994374947e-01,4.92680671114914e-17,4.58211960240755e-18,-1.68231105497676e-16,-1.11803398874990e+00,1.11022302462516e-16,-1.21136123853308e-16,2.53538066636106e-17,8.88178419700125e-16,-1.39599604810471e-17,-5.50935242191262e-17,7.25772714341823e-18,8.09016994374947e-01,8.75804297972661e-17,6.35554600906670e-17},{-4.31812442939672e-17,-1.11022302462516e-16,3.09016994374947e-01,7.95319213417573e-17,-9.91393995622428e-17,0.00000000000000e+00,-1.11803398874990e+00,-9.31654878072729e-17,9.11322387972916e-17,-0.00000000000000e+00,8.74218459219078e-16,-1.20871956352807e-16,1.62210707521434e-16,0.00000000000000e+00,8.09016994374947e-01,2.34832203558071e-17},{1.68002875787670e-16,-9.00108281395367e-17,1.31499986316849e-16,3.09016994374947e-01,-2.38434737814472e-16,-4.71549641715887e-17,-2.08260066971409e-16,-1.11803398874989e+00,1.63901328894403e-17,1.95288337460548e-16,-6.32314201647892e-17,6.96031691478364e-16,5.40417291373615e-17,-5.81225451494224e-17,1.39991500819349e-16,8.09016994374947e-01},{-5.00000000000000e-01,-1.66217541070422e-16,-1.01797970108850e-17,-3.79286206779485e-18,1.54508497187474e+00,5.63928488897144e-16,1.60634699903461e-16,-1.68742638167090e-16,-4.04508497187474e+00,-7.23573318883054e-16,-9.89869917551873e-19,2.40137976588069e-16,3.00000000000000e+00,3.25862371056332e-16,-1.49465032975025e-16,-6.76024763531834e-17},{-3.20953381081291e-16,-5.00000000000000e-01,0.00000000000000e+00,-1.23135828768772e-16,1.75993063290850e-16,1.54508497187474e+00,-8.88178419700125e-16,7.12185356409275e-16,2.70866682680208e-17,-4.04508497187474e+00,0.00000000000000e+00,-2.70866682680208e-17,-1.91751385199878e-16,3.00000000000000e+00,4.44089209850063e-16,-2.52337824650185e-16},{-4.73943045292739e-17,0.00000000000000e+00,-5.00000000000000e-01,3.19089592519251e-16,-6.78165575588495e-16,1.12895136038227e-15,1.54508497187474e+00,4.56120970663464e-16,9.47055490064095e-17,-1.40870126733237e-15,-4.04508497187474e+00,-9.47055490064095e-17,6.30854331111359e-16,-8.32667268468867e-17,3.00000000000000e+00,-6.80505014176305e-16},{-3.76857122286180e-16,9.94819828598509e-16,-1.17121716667982e-15,-5.00000000000000e-01,4.44089209850063e-16,-6.69161507634723e-16,1.39372469643533e-15,1.54508497187474e+00,0.00000000000000e+00,1.03355435391110e-16,-8.27918624191716e-16,-4.04508497187474e+00,-6.72320875638826e-17,-4.29013756354896e-16,6.05411094436203e-16,3.00000000000000e+00}};

double p_LIFT[16][16] = {{8.00000000000000e+00,-1.21014309684142e-14,-4.88498130835069e-15,0.00000000000000e+00,-2.00000000000000e+00,2.58126853225349e-15,1.24900090270330e-15,8.88178419700125e-16,9.99200722162641e-16,-1.08246744900953e-15,-6.10622663543836e-16,-2.00000000000000e+00,-4.44089209850063e-16,5.32907051820075e-15,0.00000000000000e+00,8.00000000000000e+00},
{-1.44328993201270e-15,8.00000000000000e+00,-8.88178419700125e-16,-4.44089209850063e-16,8.94427190999915e-01,-1.51267887105178e-15,-1.38777878078145e-16,5.55111512312578e-17,1.66533453693773e-16,9.29811783123569e-16,-2.00000000000000e+00,-1.31838984174237e-16,-1.11022302462516e-16,-9.71445146547012e-16,1.66533453693773e-16,-8.94427190999916e-01},{-8.88178419700125e-16,-3.71924713249427e-15,8.00000000000000e+00,2.22044604925031e-16,-8.94427190999916e-01,1.65145674912992e-15,1.66533453693773e-16,-1.66533453693773e-16,2.22044604925031e-16,-2.00000000000000e+00,2.77555756156289e-16,6.24500451351651e-17,0.00000000000000e+00,9.29811783123569e-16,-1.38777878078145e-16,8.94427190999915e-01},{-4.44089209850063e-16,4.44089209850063e-15,8.88178419700125e-16,8.00000000000000e+00,8.00000000000000e+00,-1.28785870856518e-14,-4.88498130835069e-15,0.00000000000000e+00,-2.00000000000000e+00,2.44249065417534e-15,1.55431223447522e-15,8.88178419700125e-16,9.99200722162641e-16,-1.30451205393456e-15,-6.38378239159465e-16,-2.00000000000000e+00},{-8.94427190999916e-01,1.59594559789866e-15,3.88578058618805e-16,-5.55111512312578e-17,2.91433543964104e-16,-2.00000000000000e+00,2.91433543964104e-16,5.55111512312578e-17,0.00000000000000e+00,8.74300631892311e-16,-1.80411241501588e-16,8.94427190999915e-01,-6.66133814775094e-16,-4.27435864480685e-15,8.00000000000000e+00,1.11022302462516e-16},{1.49186218934005e-16,-8.94427190999918e-01,8.32667268468867e-17,4.16333634234434e-17,-1.42247325030098e-16,8.94427190999917e-01,-2.77555756156289e-17,4.16333634234434e-17,-4.16333634234434e-17,-2.22044604925031e-16,8.94427190999916e-01,-1.04083408558608e-17,5.55111512312578e-17,2.98372437868011e-16,-8.94427190999917e-01,0.00000000000000e+00},{1.07552855510562e-16,3.74700270810990e-16,-8.94427190999917e-01,-1.38777878078145e-17,1.66533453693773e-16,-8.94427190999917e-01,-6.93889390390723e-18,-1.38777878078145e-17,-1.80411241501588e-16,8.94427190999917e-01,6.93889390390723e-17,2.08166817117217e-17,-2.77555756156289e-17,-3.05311331771918e-16,8.94427190999916e-01,2.42861286636753e-17},{-2.77555756156289e-17,-7.21644966006352e-16,1.11022302462516e-16,-8.94427190999916e-01,-1.38777878078145e-15,8.00000000000000e+00,-1.72084568816899e-15,-2.22044604925031e-16,8.94427190999916e-01,-1.66533453693773e-15,-1.66533453693773e-16,5.55111512312578e-17,1.66533453693773e-16,9.29811783123569e-16,-2.00000000000000e+00,-8.32667268468867e-17},{8.94427190999916e-01,-1.66533453693773e-15,-1.52655665885959e-16,5.55111512312578e-17,2.01227923213310e-16,9.15933995315754e-16,-2.00000000000000e+00,-5.55111512312578e-17,-5.55111512312578e-17,-7.49400541621981e-16,1.11022302462516e-16,-8.94427190999916e-01,-1.33226762955019e-15,8.00000000000001e+00,-1.77635683940025e-15,-2.22044604925031e-16},{-1.87350135405495e-16,8.94427190999917e-01,5.55111512312578e-17,1.38777878078145e-17,-4.51028103753970e-17,-3.12250225675825e-16,8.94427190999916e-01,2.77555756156289e-17,1.11022302462516e-16,3.60822483003176e-16,-8.94427190999917e-01,-1.73472347597681e-17,1.66533453693773e-16,-8.94427190999917e-01,6.93889390390723e-18,-1.73472347597681e-17},{-4.85722573273506e-17,-2.15105711021124e-16,8.94427190999916e-01,-1.38777878078145e-17,5.20417042793042e-17,2.91433543964104e-16,-8.94427190999917e-01,0.00000000000000e+00,1.52655665885959e-16,-8.94427190999918e-01,9.02056207507940e-17,3.81639164714898e-17,-1.38777878078145e-16,8.94427190999917e-01,-2.08166817117217e-17,3.81639164714898e-17},{0.00000000000000e+00,8.88178419700125e-16,-1.66533453693773e-16,8.94427190999915e-01,-7.77156117237610e-16,-4.16333634234434e-15,8.00000000000000e+00,2.22044604925031e-16,-8.94427190999916e-01,1.60982338570648e-15,3.33066907387547e-16,-5.55111512312578e-17,2.77555756156289e-16,-2.00000000000000e+00,2.63677968348475e-16,6.24500451351651e-17},{-2.00000000000000e+00,2.38697950294409e-15,1.47104550762833e-15,8.88178419700125e-16,9.99200722162641e-16,-1.22124532708767e-15,-6.66133814775094e-16,-2.00000000000000e+00,-4.44089209850063e-16,4.44089209850063e-15,7.77156117237610e-16,8.00000000000000e+00,8.00000000000000e+00,-1.28785870856518e-14,-5.32907051820075e-15,-2.22044604925031e-16},{2.42861286636753e-16,-2.00000000000000e+00,3.05311331771918e-16,5.55111512312578e-17,0.00000000000000e+00,9.43689570931383e-16,-1.11022302462516e-16,8.94427190999915e-01,-8.88178419700125e-16,-3.88578058618805e-15,8.00000000000000e+00,8.32667268468867e-17,-8.94427190999916e-01,1.66533453693773e-15,1.66533453693773e-16,-1.11022302462516e-16},{2.01227923213310e-16,9.29811783123569e-16,-2.00000000000000e+00,-1.11022302462516e-16,-8.32667268468867e-17,-9.43689570931383e-16,1.66533453693773e-16,-8.94427190999916e-01,-1.55431223447522e-15,8.00000000000001e+00,-9.43689570931383e-16,-3.33066907387547e-16,8.94427190999915e-01,-1.49880108324396e-15,-1.66533453693773e-16,5.55111512312578e-17},{9.99200722162641e-16,-1.11022302462516e-15,-6.66133814775094e-16,-2.00000000000000e+00,-4.44089209850063e-16,5.32907051820075e-15,0.00000000000000e+00,8.00000000000000e+00,8.00000000000000e+00,-1.19904086659517e-14,-5.32907051820075e-15,-2.22044604925031e-16,-2.00000000000000e+00,2.55351295663786e-15,1.44328993201270e-15,8.32667268468867e-16}};

int p_Fmask[4][4] = {{0,1,2,3},
{3,7,11,15},{15,14,13,12},{12,8,4,0}};

#endif
