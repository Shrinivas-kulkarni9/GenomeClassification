function [feature_mtx,labels,num_labels] = readData(path_to_data,encoding_scheme)
    %READDATA reads DNA sequences specified by path_to_data input, maps
    %them to numbers using encoding_scheme provided as input, median
    %normalises them, and returns output as described below
    %   Inputs:
    %       path_to_data: relative path of dataset to read. To be given as
    %           a string
    %       encoding_scheme: encoding scheme to be used to get numerical
    %           sequence corresponding to input DNA_seq.
    %           Possible values:
    %             1)PP: Purine/Pyrimidine, [T,C,A,G] -> [1,1,-1,-1]
    %             2)Just-A: Just-A, [T,C,A,G] -> [0,0,1,0]
    %             3)Real: Real, [T,C,A,G] -> [-1.5,0.5,1.5,-0.5]
    %   Outputs:
    %       feature_mtx: each column contains a DNA sequence, each of which
    %           is median normalised
    %       labels: a vector whose ith element is the label of DNA sequence
    %           in ith column of feature_mtx
    %       num_labels: number of distinct labels in labels. Hence, labels
    %           contains values from 1 to num_labels (both included)

    % read all individual sequences -> make every DNA_seq numerical (while
    % reading) -> calc median -> normalise each sequence -> make matrix and
    % return it
    direc = dir(path_to_data);
    num_labels = size(direc,1)-3;
    num_sequences = 0;
    labels = 0;
    lengths = 0;
    labels_cell_array = cell(num_labels,1);
    for i=1:num_labels
        direc_this_label = dir(strcat(path_to_data,direc(i+3).name,"/*.txt"));
        num_seq_this_label = size(direc_this_label,1);
        
        labels(num_sequences+1:num_sequences+num_seq_this_label) = i;
        
        seq_cell_array = cell(num_seq_this_label,1);
        for j=1:num_seq_this_label
            file_path = strcat(path_to_data,direc(i+3).name,"/",direc_this_label(j).name);
            fileID = fopen(file_path,'r');
            [file_content,len_this_seq] = fscanf(fileID,'%c');
            new_line_strt_indices = regexp(file_content,'\n');
            file_content(new_line_strt_indices) = [];
            len_this_seq = len_this_seq - size(new_line_strt_indices,2);
            file_content = file_content(new_line_strt_indices(1):len_this_seq);
            len_this_seq = size(file_content,2);
            
            seq_cell_array{j} = file_content;
            lengths(num_sequences+j) = len_this_seq;
        end
        labels_cell_array{i} = seq_cell_array;
        num_sequences = num_sequences + num_seq_this_label;
    end
    
    median_length = floor(median(lengths));
    feature_mtx = zeros([median_length,num_sequences]);
    
    k = 0; % used only for some temp indexing, goes from 1 to num_sequences
    for i = 1:num_labels
        this_label_cell = labels_cell_array{i};
        for j = 1:size(this_label_cell,1)
            this_seq_num = DNA_SeqToNum(this_label_cell{j},encoding_scheme);
            if size(this_seq_num,2) > median_length
                % truncate
                this_seq_num = this_seq_num(1:median_length);
            elseif size(this_seq_num,2) == median_length
                % do nothing
            else
                % up-sample
                this_seq_mirror = -flip(this_seq_num);
                while(size(this_seq_num,2) < median_length)
                    if(size(this_seq_num,2)+size(this_seq_mirror,2) <= median_length)
                        this_seq_num(size(this_seq_num,2)+1:size(this_seq_num,2)+size(this_seq_mirror,2))...
                            = this_seq_mirror;
                    else
                        this_seq_num(size(this_seq_num,2)+1:median_length) = this_seq_mirror...
                            (1:median_length-size(this_seq_num,2));
                    end
                    this_seq_mirror = -flip(this_seq_mirror);
                end
            end
            feature_mtx(:,k+1) = this_seq_num;
            k = k+1;
        end
    end

end

