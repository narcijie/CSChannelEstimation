function [output,count,pilot_sequence]=Insert_pilot(pilot_inter,pilot_symbol_bit,map_out_block,method) 

switch method
    case 'QAM16'
        pilot_symbol=QAM16(pilot_symbol_bit);     
    case 'QPSK'
        s=(pilot_symbol_bit.*2-1)/sqrt(2);
        pilot_symbol=complex(s(1),s(2));
end

        [N,NL]=size(map_out_block);
        output=zeros(N,(NL+fix(NL/pilot_inter)));
        pilot_sequence=pilot_symbol*ones(N,1);
        
        count=0;
        i=1;
        while i<(NL+fix(NL/pilot_inter))
            output(:,i)=pilot_sequence;
            count=count+1;
            if count*pilot_inter<=NL
                output(:,(i+1):(i+pilot_inter))=map_out_block(:,((count-1)*pilot_inter+1):count*pilot_inter);
            else
                output(:,(i+1):(i+pilot_inter+NL-count*pilot_inter))=map_out_block(:,((count-1)*pilot_inter+1):NL);
            end
            i=i+pilot_inter+1;
        end
end

    