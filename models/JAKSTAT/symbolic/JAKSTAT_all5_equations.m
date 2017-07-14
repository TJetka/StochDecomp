function R = f(t,y,p)

R = [
     [ (2*p(4)*y(13) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*JAKSTAT_stimulus(t)) ];
     [ (exp(y(17) - p(8)/p(9))*p(1)*y(1)*JAKSTAT_stimulus(t) - 2*p(2)*y(2)^2) ];
     [ (p(2)*y(2)^2 - p(3)*y(3)) ];
     [ (p(3)*y(3) - p(4)*y(4)) ];
     [ (p(4)*y(4) - p(4)*y(5)) ];
     [ (p(4)*y(5) - p(4)*y(6)) ];
     [ (p(4)*y(6) - p(4)*y(7)) ];
     [ (p(4)*y(7) - p(4)*y(8)) ];
     [ (p(4)*y(8) - p(4)*y(9)) ];
     [ (p(4)*y(9) - p(4)*y(10)) ];
     [ (p(4)*y(10) - p(4)*y(11)) ];
     [ (p(4)*y(11) - p(4)*y(12)) ];
     [ (p(4)*y(12) - p(4)*y(13)) ];
     [ (exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*JAKSTAT_stimulus(t) - 2*p(3)*p(6)*y(3)) ];
     [ (2*p(4)*p(7)*y(13) - 2*p(3)*p(7)*y(3)) ];
     [ (p(3)*y(3) - p(4)*y(13)) ];
     [ (p(8) - p(9)*y(17)) ];
     [ (4*p(4)*y(46) - 2*exp(y(17) - p(8)/p(9))*p(1)*y(18)*JAKSTAT_stimulus(t) - 2*exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(50)*JAKSTAT_stimulus(t)) ];
     [ (2*exp(y(17) - p(8)/p(9))*p(1)*y(35)*JAKSTAT_stimulus(t) - 8*p(2)*y(2)*y(19) + 2*exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(65)*JAKSTAT_stimulus(t)) ];
     [ (4*p(2)*y(2)*y(51) - 2*p(3)*y(20)) ];
     [ (2*p(3)*y(66) - 2*p(4)*y(21)) ];
     [ (p(4)*y(5) - 2*p(4)*y(22) + 2*p(4)*y(80)) ];
     [ (p(4)*y(5) - 2*p(4)*y(23) + 2*p(4)*y(93)) ];
     [ (2*p(4)*y(105) - 2*p(4)*y(24)) ];
     [ (2*p(4)*y(116) - 2*p(4)*y(25)) ];
     [ (2*p(4)*y(126) - 2*p(4)*y(26)) ];
     [ (2*p(4)*y(135) - 2*p(4)*y(27)) ];
     [ (2*p(4)*y(143) - 2*p(4)*y(28)) ];
     [ (2*p(4)*y(150) - 2*p(4)*y(29)) ];
     [ (2*p(4)*y(156) - 2*p(4)*y(30)) ];
     [ (2*exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(47)*JAKSTAT_stimulus(t) - 4*p(3)*p(6)*y(76) + 2*exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(167)*JAKSTAT_stimulus(t)) ];
     [ (4*p(4)*p(7)*y(162) - 4*p(3)*p(7)*y(77)) ];
     [ (2*p(3)*y(78) - 2*p(4)*y(163)) ];
     [ (-2*p(9)*y(34)) ];
     [ (2*p(4)*y(61) - 4*p(2)*y(2)*y(35) + exp(y(17) - p(8)/p(9))*p(1)*y(18)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(35)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(50)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(65)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*y(75) - p(3)*y(36) + 2*p(2)*y(2)*y(35) - exp(y(17) - p(8)/p(9))*p(1)*y(36)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(79)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(36) - p(4)*y(37) + 2*p(4)*y(88) - exp(y(17) - p(8)/p(9))*p(1)*y(37)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(92)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(37) - p(4)*y(38) + 2*p(4)*y(100) - exp(y(17) - p(8)/p(9))*p(1)*y(38)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(104)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(38) - p(4)*y(39) + 2*p(4)*y(111) - exp(y(17) - p(8)/p(9))*p(1)*y(39)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(115)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(39) - p(4)*y(40) + 2*p(4)*y(121) - exp(y(17) - p(8)/p(9))*p(1)*y(40)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(125)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(40) - p(4)*y(41) + 2*p(4)*y(130) - exp(y(17) - p(8)/p(9))*p(1)*y(41)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(134)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(41) - p(4)*y(42) + 2*p(4)*y(138) - exp(y(17) - p(8)/p(9))*p(1)*y(42)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(142)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(42) - p(4)*y(43) + 2*p(4)*y(145) - exp(y(17) - p(8)/p(9))*p(1)*y(43)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(149)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(43) - p(4)*y(44) + 2*p(4)*y(151) - exp(y(17) - p(8)/p(9))*p(1)*y(44)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(155)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(44) - p(4)*y(45) + 2*p(4)*y(156) - exp(y(17) - p(8)/p(9))*p(1)*y(45)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(160)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*y(30) + p(4)*y(45) - p(4)*y(46) - exp(y(17) - p(8)/p(9))*p(1)*y(46)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(164)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*y(161) - 2*p(3)*p(6)*y(36) - exp(y(17) - p(8)/p(9))*p(1)*y(47)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(18)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(167)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(50)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*y(162) - 2*p(3)*p(7)*y(36) + 2*p(4)*p(7)*y(46) - exp(y(17) - p(8)/p(9))*p(1)*y(48)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(169)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(36) - p(4)*y(46) + 2*p(4)*y(163) - exp(y(17) - p(8)/p(9))*p(1)*y(49)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(170)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*y(164) - p(9)*y(50) - exp(y(17) - p(8)/p(9))*p(1)*y(50)*JAKSTAT_stimulus(t) - exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(34)*JAKSTAT_stimulus(t)) ];
     [ (2*p(2)*y(2)*y(19) - p(3)*y(51) - 4*p(2)*y(2)*y(51) + exp(y(17) - p(8)/p(9))*p(1)*y(36)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(79)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(51) - p(4)*y(52) - 4*p(2)*y(2)*y(52) + exp(y(17) - p(8)/p(9))*p(1)*y(37)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(92)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(52) - p(4)*y(53) - 4*p(2)*y(2)*y(53) + exp(y(17) - p(8)/p(9))*p(1)*y(38)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(104)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(53) - p(4)*y(54) - 4*p(2)*y(2)*y(54) + exp(y(17) - p(8)/p(9))*p(1)*y(39)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(115)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(54) - p(4)*y(55) - 4*p(2)*y(2)*y(55) + exp(y(17) - p(8)/p(9))*p(1)*y(40)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(125)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(55) - p(4)*y(56) - 4*p(2)*y(2)*y(56) + exp(y(17) - p(8)/p(9))*p(1)*y(41)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(134)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(56) - p(4)*y(57) - 4*p(2)*y(2)*y(57) + exp(y(17) - p(8)/p(9))*p(1)*y(42)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(142)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(57) - p(4)*y(58) - 4*p(2)*y(2)*y(58) + exp(y(17) - p(8)/p(9))*p(1)*y(43)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(149)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(58) - p(4)*y(59) - 4*p(2)*y(2)*y(59) + exp(y(17) - p(8)/p(9))*p(1)*y(44)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(155)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(59) - p(4)*y(60) - 4*p(2)*y(2)*y(60) + exp(y(17) - p(8)/p(9))*p(1)*y(45)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(160)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(60) - p(4)*y(61) - 4*p(2)*y(2)*y(61) + exp(y(17) - p(8)/p(9))*p(1)*y(46)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(164)*JAKSTAT_stimulus(t)) ];
     [ (exp(y(17) - p(8)/p(9))*p(1)*y(47)*JAKSTAT_stimulus(t) - 4*p(2)*y(2)*y(62) - 2*p(3)*p(6)*y(51) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(35)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(167)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(65)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*p(7)*y(61) - 2*p(3)*p(7)*y(51) - 4*p(2)*y(2)*y(63) + exp(y(17) - p(8)/p(9))*p(1)*y(48)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(169)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(51) - p(4)*y(61) - 4*p(2)*y(2)*y(64) + exp(y(17) - p(8)/p(9))*p(1)*y(49)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(170)*JAKSTAT_stimulus(t)) ];
     [ (exp(y(17) - p(8)/p(9))*p(1)*y(50)*JAKSTAT_stimulus(t) - 4*p(2)*y(2)*y(65) - p(9)*y(65) + exp(y(17) - p(8)/p(9))*p(1)*y(1)*y(34)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(20) - p(3)*y(66) - p(4)*y(66) + 2*p(2)*y(2)*y(52)) ];
     [ (p(4)*y(66) - p(3)*y(67) - p(4)*y(67) + 2*p(2)*y(2)*y(53)) ];
     [ (p(4)*y(67) - p(3)*y(68) - p(4)*y(68) + 2*p(2)*y(2)*y(54)) ];
     [ (p(4)*y(68) - p(3)*y(69) - p(4)*y(69) + 2*p(2)*y(2)*y(55)) ];
     [ (p(4)*y(69) - p(3)*y(70) - p(4)*y(70) + 2*p(2)*y(2)*y(56)) ];
     [ (p(4)*y(70) - p(3)*y(71) - p(4)*y(71) + 2*p(2)*y(2)*y(57)) ];
     [ (p(4)*y(71) - p(3)*y(72) - p(4)*y(72) + 2*p(2)*y(2)*y(58)) ];
     [ (p(4)*y(72) - p(3)*y(73) - p(4)*y(73) + 2*p(2)*y(2)*y(59)) ];
     [ (p(4)*y(73) - p(3)*y(74) - p(4)*y(74) + 2*p(2)*y(2)*y(60)) ];
     [ (p(4)*y(74) - p(3)*y(75) - p(4)*y(75) + 2*p(2)*y(2)*y(61)) ];
     [ (2*p(2)*y(2)*y(62) - 2*p(3)*p(6)*y(20) - p(3)*y(76) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(36)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(79)*JAKSTAT_stimulus(t)) ];
     [ (2*p(4)*p(7)*y(75) - 2*p(3)*p(7)*y(20) - p(3)*y(77) + 2*p(2)*y(2)*y(63)) ];
     [ (p(3)*y(20) - p(4)*y(75) - p(3)*y(78) + 2*p(2)*y(2)*y(64)) ];
     [ (2*p(2)*y(2)*y(65) - p(9)*y(79) - p(3)*y(79)) ];
     [ (p(4)*y(21) + p(3)*y(67) - 2*p(4)*y(80)) ];
     [ (p(3)*y(68) + p(4)*y(80) - 2*p(4)*y(81)) ];
     [ (p(3)*y(69) + p(4)*y(81) - 2*p(4)*y(82)) ];
     [ (p(3)*y(70) + p(4)*y(82) - 2*p(4)*y(83)) ];
     [ (p(3)*y(71) + p(4)*y(83) - 2*p(4)*y(84)) ];
     [ (p(3)*y(72) + p(4)*y(84) - 2*p(4)*y(85)) ];
     [ (p(3)*y(73) + p(4)*y(85) - 2*p(4)*y(86)) ];
     [ (p(3)*y(74) + p(4)*y(86) - 2*p(4)*y(87)) ];
     [ (p(3)*y(75) + p(4)*y(87) - 2*p(4)*y(88)) ];
     [ (p(3)*y(76) - p(4)*y(89) - 2*p(3)*p(6)*y(66) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(37)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(92)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(77) - p(4)*y(90) - 2*p(3)*p(7)*y(66) + 2*p(4)*p(7)*y(88)) ];
     [ (p(3)*y(66) + p(3)*y(78) - p(4)*y(88) - p(4)*y(91)) ];
     [ (p(3)*y(79) - p(4)*y(92) - p(9)*y(92)) ];
     [ (p(4)*y(22) - p(4)*y(5) + p(4)*y(81) - 2*p(4)*y(93)) ];
     [ (p(4)*y(82) + p(4)*y(93) - 2*p(4)*y(94)) ];
     [ (p(4)*y(83) + p(4)*y(94) - 2*p(4)*y(95)) ];
     [ (p(4)*y(84) + p(4)*y(95) - 2*p(4)*y(96)) ];
     [ (p(4)*y(85) + p(4)*y(96) - 2*p(4)*y(97)) ];
     [ (p(4)*y(86) + p(4)*y(97) - 2*p(4)*y(98)) ];
     [ (p(4)*y(87) + p(4)*y(98) - 2*p(4)*y(99)) ];
     [ (p(4)*y(88) + p(4)*y(99) - 2*p(4)*y(100)) ];
     [ (p(4)*y(89) - p(4)*y(101) - 2*p(3)*p(6)*y(67) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(38)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(104)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(90) - p(4)*y(102) - 2*p(3)*p(7)*y(67) + 2*p(4)*p(7)*y(100)) ];
     [ (p(3)*y(67) + p(4)*y(91) - p(4)*y(100) - p(4)*y(103)) ];
     [ (p(4)*y(92) - p(4)*y(104) - p(9)*y(104)) ];
     [ (p(4)*y(23) + p(4)*y(94) - 2*p(4)*y(105)) ];
     [ (p(4)*y(95) + p(4)*y(105) - 2*p(4)*y(106)) ];
     [ (p(4)*y(96) + p(4)*y(106) - 2*p(4)*y(107)) ];
     [ (p(4)*y(97) + p(4)*y(107) - 2*p(4)*y(108)) ];
     [ (p(4)*y(98) + p(4)*y(108) - 2*p(4)*y(109)) ];
     [ (p(4)*y(99) + p(4)*y(109) - 2*p(4)*y(110)) ];
     [ (p(4)*y(100) + p(4)*y(110) - 2*p(4)*y(111)) ];
     [ (p(4)*y(101) - p(4)*y(112) - 2*p(3)*p(6)*y(68) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(39)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(115)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(102) - p(4)*y(113) - 2*p(3)*p(7)*y(68) + 2*p(4)*p(7)*y(111)) ];
     [ (p(3)*y(68) + p(4)*y(103) - p(4)*y(111) - p(4)*y(114)) ];
     [ (p(4)*y(104) - p(4)*y(115) - p(9)*y(115)) ];
     [ (p(4)*y(24) + p(4)*y(106) - 2*p(4)*y(116)) ];
     [ (p(4)*y(107) + p(4)*y(116) - 2*p(4)*y(117)) ];
     [ (p(4)*y(108) + p(4)*y(117) - 2*p(4)*y(118)) ];
     [ (p(4)*y(109) + p(4)*y(118) - 2*p(4)*y(119)) ];
     [ (p(4)*y(110) + p(4)*y(119) - 2*p(4)*y(120)) ];
     [ (p(4)*y(111) + p(4)*y(120) - 2*p(4)*y(121)) ];
     [ (p(4)*y(112) - p(4)*y(122) - 2*p(3)*p(6)*y(69) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(40)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(125)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(113) - p(4)*y(123) - 2*p(3)*p(7)*y(69) + 2*p(4)*p(7)*y(121)) ];
     [ (p(3)*y(69) + p(4)*y(114) - p(4)*y(121) - p(4)*y(124)) ];
     [ (p(4)*y(115) - p(4)*y(125) - p(9)*y(125)) ];
     [ (p(4)*y(25) + p(4)*y(117) - 2*p(4)*y(126)) ];
     [ (p(4)*y(118) + p(4)*y(126) - 2*p(4)*y(127)) ];
     [ (p(4)*y(119) + p(4)*y(127) - 2*p(4)*y(128)) ];
     [ (p(4)*y(120) + p(4)*y(128) - 2*p(4)*y(129)) ];
     [ (p(4)*y(121) + p(4)*y(129) - 2*p(4)*y(130)) ];
     [ (p(4)*y(122) - p(4)*y(131) - 2*p(3)*p(6)*y(70) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(41)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(134)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(123) - p(4)*y(132) - 2*p(3)*p(7)*y(70) + 2*p(4)*p(7)*y(130)) ];
     [ (p(3)*y(70) + p(4)*y(124) - p(4)*y(130) - p(4)*y(133)) ];
     [ (p(4)*y(125) - p(4)*y(134) - p(9)*y(134)) ];
     [ (p(4)*y(26) + p(4)*y(127) - 2*p(4)*y(135)) ];
     [ (p(4)*y(128) + p(4)*y(135) - 2*p(4)*y(136)) ];
     [ (p(4)*y(129) + p(4)*y(136) - 2*p(4)*y(137)) ];
     [ (p(4)*y(130) + p(4)*y(137) - 2*p(4)*y(138)) ];
     [ (p(4)*y(131) - p(4)*y(139) - 2*p(3)*p(6)*y(71) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(42)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(142)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(132) - p(4)*y(140) - 2*p(3)*p(7)*y(71) + 2*p(4)*p(7)*y(138)) ];
     [ (p(3)*y(71) + p(4)*y(133) - p(4)*y(138) - p(4)*y(141)) ];
     [ (p(4)*y(134) - p(4)*y(142) - p(9)*y(142)) ];
     [ (p(4)*y(27) + p(4)*y(136) - 2*p(4)*y(143)) ];
     [ (p(4)*y(137) + p(4)*y(143) - 2*p(4)*y(144)) ];
     [ (p(4)*y(138) + p(4)*y(144) - 2*p(4)*y(145)) ];
     [ (p(4)*y(139) - p(4)*y(146) - 2*p(3)*p(6)*y(72) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(43)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(149)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(140) - p(4)*y(147) - 2*p(3)*p(7)*y(72) + 2*p(4)*p(7)*y(145)) ];
     [ (p(3)*y(72) + p(4)*y(141) - p(4)*y(145) - p(4)*y(148)) ];
     [ (p(4)*y(142) - p(4)*y(149) - p(9)*y(149)) ];
     [ (p(4)*y(28) + p(4)*y(144) - 2*p(4)*y(150)) ];
     [ (p(4)*y(145) + p(4)*y(150) - 2*p(4)*y(151)) ];
     [ (p(4)*y(146) - p(4)*y(152) - 2*p(3)*p(6)*y(73) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(44)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(155)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(147) - p(4)*y(153) - 2*p(3)*p(7)*y(73) + 2*p(4)*p(7)*y(151)) ];
     [ (p(3)*y(73) + p(4)*y(148) - p(4)*y(151) - p(4)*y(154)) ];
     [ (p(4)*y(149) - p(4)*y(155) - p(9)*y(155)) ];
     [ (p(4)*y(29) + p(4)*y(151) - 2*p(4)*y(156)) ];
     [ (p(4)*y(152) - p(4)*y(157) - 2*p(3)*p(6)*y(74) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(45)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(160)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(153) - p(4)*y(158) - 2*p(3)*p(7)*y(74) + 2*p(4)*p(7)*y(156)) ];
     [ (p(3)*y(74) + p(4)*y(154) - p(4)*y(156) - p(4)*y(159)) ];
     [ (p(4)*y(155) - p(4)*y(160) - p(9)*y(160)) ];
     [ (p(4)*y(157) - p(4)*y(161) - 2*p(3)*p(6)*y(75) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(46)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(164)*JAKSTAT_stimulus(t)) ];
     [ (p(4)*y(158) - p(4)*y(162) + 2*p(4)*p(7)*y(30) - 2*p(3)*p(7)*y(75)) ];
     [ (p(3)*y(75) - p(4)*y(30) + p(4)*y(159) - p(4)*y(163)) ];
     [ (p(4)*y(160) - p(4)*y(164) - p(9)*y(164)) ];
     [ (2*p(4)*p(7)*y(161) - 2*p(3)*p(7)*y(76) - 2*p(3)*p(6)*y(77) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(48)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(169)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(76) - p(4)*y(161) - 2*p(3)*p(6)*y(78) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(49)*JAKSTAT_stimulus(t) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(170)*JAKSTAT_stimulus(t)) ];
     [ (exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(50)*JAKSTAT_stimulus(t) - 2*p(3)*p(6)*y(79) - p(9)*y(167) + exp(y(17) - p(8)/p(9))*p(1)*p(6)*y(1)*y(34)*JAKSTAT_stimulus(t)) ];
     [ (p(3)*y(77) - p(4)*y(162) - 2*p(3)*p(7)*y(78) + 2*p(4)*p(7)*y(163)) ];
     [ (2*p(4)*p(7)*y(164) - 2*p(3)*p(7)*y(79) - p(9)*y(169)) ];
     [ (p(3)*y(79) - p(4)*y(164) - p(9)*y(170)) ];
];
