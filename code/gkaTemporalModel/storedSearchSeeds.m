function pMRI0 = storedSearchSeeds(whichSub,modelType)

% Stored p0 vectors
switch whichSub
    case 1
        switch modelType
            case 'stimulus' % fVal = 6.04
                pMRI0 = [ ...
                    0.8591833733, 29.8358988389, 44.6459064931,  ...
                    0.9417344220, 0.1241009757, 0.6671111792, 5.9140191451, 0.5713267326, 1.5196739733,  ...
                    29.7649686038, 0.9994013488, 0.7555992424, 0.7549213707, 0.7440302610, 0.6354759276, 0.6304316044, 0.0951708704, 0.1724046767, 0.5084658861, 0.4943162948, 0.5664237291, 1.5087736771,  ...
                    4.0573334694, 0.9996716380, 0.5949135602, 0.4713958919, 0.4594216347, 0.3954613566, 0.2865379333, 7.2781368718, 21.3026127368, 28.1745418012, 18.1641720384, 4.8454543650, 2.8494692817,  ...
                    5.9936985373, 0.7544339120, 0.7502034307, 0.7396474838, 0.5155055821, -0.2870200813, -0.9570114553, 0.6806548536, 1.1499288455, 1.9817052633, 1.7818458229, 0.7818616405, 0.9763333276,  ...
                    ];
            case 'cell' % fVal = 7.28
                pMRI0 = [ ...
                    0.8706967235, 46.8667360544, 45.3693057299,  ...
                    0.8982335091, 0.1235109568, 0.5981354117, 5.7011002302, 0.5097818732, 2.1928362846,  ...
                    23.3449769020, 0.9256205559, 0.6311378479, 0.5759220123, 0.3442761421, 0.6047307014, 0.3031565666, 0.0979398489, 0.1821182966, 0.5207353830, 0.4823014736, 0.3225244284, 0.6256033182,  ...
                    3.0001521111, 0.8546582222, 0.6589594841, 0.4955028534, 0.5030731201, 0.3181489944, 0.5156174660, 8.1160496473, 24.4760655165, 27.6791831255, 18.7499837875, 3.9993255138, 3.9414112568,  ...
                    7.8980207443, 0.6555871010, 0.7443149567, 0.5826660156, 0.3865708351, -0.7406884193, -0.9994610786, 1.5763520002, 2.1527512074, 2.1747303009, 2.1960763931, 0.8695276976, 1.1486126184,  ...
                    ];
        end
    case 2
        switch modelType
            case 'stimulus' % fVal = 4.73
                pMRI0 = [ ...
                    0.7699677318, 29.3925114274, 31.1932458878,  ...
                    0.9933977216, 0.3457751572, 0.1686610967, 1.0920037925, 0.9703282386, 1.3071486056,  ...
                    29.9977755547, 0.4511427402, 0.4268976688, 0.2643770933, -0.0733599424, -0.1393993616, -0.1400948286, 0.1261054575, 0.2925901711, 0.9176430404, 1.0880668461, 2.8610488474, 11.0773752332,  ...
                    3.5546475649, 0.4622892141, 0.2796707392, 0.0897983074, 0.0518857241, -0.0445317507, -0.1942487717, 3.6534522474, 15.5568882525, 25.3146409988, 21.0678004920, 14.5492617190, 9.3953396976,  ...
                    12.2112888098, 0.9790496826, 0.9430859804, 0.8439409018, 0.5629886866, 0.1973487616, -0.7807117701, 0.5094105601, 1.1186485887, 2.1984906495, 2.3175396919, 2.7600489259, 3.7193806767,  ...
                    ];
            case 'cell' % fVal = 5.05
                pMRI0 = [ ...
                    0.7550306804, 32.1684170812, 31.6757655814,  ...
                    0.9991397411, 0.3228568882, 0.1432418399, 1.2317744419, 0.9998323098, 2.4885661900,  ...
                    29.9996089935, 0.4785507083, 0.4098671079, 0.2300033510, -0.1156228304, -0.0637271583, 0.0754283726, 0.1274390519, 0.2922273353, 0.9428859130, 1.1875284538, 2.9240133017, 10.5693342686,  ...
                    3.5871683061, 0.5028908908, 0.3609825552, 0.0508957744, 0.0798755229, 0.0334557295, -0.4830626369, 5.8000357226, 16.3618021458, 24.1258842200, 22.0267170370, 15.5494596064, 8.0120950043,  ...
                    12.9364675283, 0.9957849979, 0.9952014923, 0.9995374918, 0.6622459412, 0.2853996217, -0.7530252457, 1.8948657289, 2.9901290163, 4.0457154587, 3.3188655823, 2.4391760603, 3.0404706597,  ...
                    ];
        end


end