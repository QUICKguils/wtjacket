classdef HarmonicLoad < Load
	% HARMONICLOAD  Represent a harmonic structural load.

	properties (SetAccess = protected)
		amplitude      (1,  1) double {mustBeReal}  % Amplitude [N].
		frequencyHertz (1,  1) double {mustBeReal}  % Frequency [Hz].
		frequencyRad   (1,  1) double {mustBeReal}  % Frequency [rad/s].
	end

	methods
		function hLoad = HarmonicLoad(node, direction, amplitude, frequencyHertz)
			timeEvolution = @(t) amplitude*sin(frequencyHertz*2*pi*t);

			hLoad = hLoad@Load(node, direction, timeEvolution);

			hLoad.amplitude      = amplitude;
			hLoad.frequencyHertz = frequencyHertz;
			hLoad.frequencyRad   = frequencyHertz * 2*pi;
		end
	end
end
