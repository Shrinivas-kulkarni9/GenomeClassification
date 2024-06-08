function [num_seq] = DNA_SeqToNum(DNA_seq,encoding_scheme)
    %DNA_SEQTONUM Takes as input a DNA sequence and returns the encoded
    %numerical vector using input encoding scheme
    %   Inputs:
    %       DNA_seq: DNA sequence input as a char array containing A,C,G,T
    %       encoding_scheme: encoding scheme to be used to get numerical
    %           sequence corresponding to input DNA_seq.
    %           Possible values:
    %             1)PP: Purine/Pyrimidine, [T,C,A,G] -> [1,1,-1,-1]
    %             2)Just-A: Just-A, [T,C,A,G] -> [0,0,1,0]
    %             3)Real: Real, [T,C,A,G] -> [-1.5,0.5,1.5,-0.5]
    %   Outputs:
    %       num_seq: encoded DNA_seq using encoding_scheme
    num_seq = zeros(size(DNA_seq));
    switch encoding_scheme
        case "PP"
            num_seq(DNA_seq == 'T') = 1;
            num_seq(DNA_seq == 'C') = 1;
            num_seq(DNA_seq == 'A') = -1;
            num_seq(DNA_seq == 'G') = -1;
        case "Just-A"
            num_seq(DNA_seq == 'A') = 1;
        case "Real"
            num_seq(DNA_seq == 'T') = -1.5;
            num_seq(DNA_seq == 'C') = 0.5;
            num_seq(DNA_seq == 'A') = 1.5;
            num_seq(DNA_seq == 'G') = -0.5;
        otherwise
            disp("Please enter one of encoding schemes given in function description.");
    end
end

