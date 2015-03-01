%% Compile mex function Doptim

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
fprintf('compilation of Doptim mex function (octave %d)\n', isOctave);

% define sources
lst={'src/anyoption.cpp', 'src/Deff.cpp', 'src/oaoptions.cpp', 'src/tools.cpp', 'src/arrayproperties.cpp', 'src/mathtools.cpp', 'src/pareto.cpp', 'src/InfInt.cpp','src/strength.cpp', 'src/md5.cpp', 'src/arraytools.cpp'};
lst =cat(1,lst(:), 'bitarray/bit_array.cpp');

% add lmc code
lst =cat(1,lst(:), 'src/lmc.cpp', 'src/extend.cpp', 'src/nonroot.cpp');

fprintf('  compiling Doptim.cpp with %d source files\n', length(lst));

args={'-Isrc', '-I.'};

if isOctave
    args{end+1}='-DNOOPENMP';
else
    %args{end+1}='-DDOOPENMP';
end
args{end+1}='-v';

% compile
mex('mex/Doptim.cpp', '-o', 'mex/Doptim', lst{:}, args{:},  '-DNOZLIB'  )
fprintf('compilation of Doptim complete\n');

return;

% test: mkoctfile --mex mex/Doptim.cpp -I. -Isrc src/*cpp -o mex/Doptim

%% Test the resulting mex file

[A, d] = Doptim(40, [2,2,2,2,2,2], 10, [1,2,0],1);


