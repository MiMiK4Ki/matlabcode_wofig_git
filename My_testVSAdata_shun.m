function [iqdata_vsa1, iqdata_vsa2] = My_testVSAdata_shun(vsa, Ns, tracenum)
    % Close and reopen VSA connection
    fclose(vsa);
    vsa.TimeOut = 10;
    vsa.byteOrder = 'bigEndian';
    vsa.inputbuffersize = 50e6;
    vsa.outputbuffersize = 50e6;
    fopen(vsa);
    % Initialize variables
    DataFormat = 'real64';
    % Ensure tracenum contains two trace numbers
% % %     if numel(tracenum) ~= 2
% % %         error('tracenum must contain exactly two trace numbers');
% % %     end
    trace1 = tracenum(1);
% % % % %     trace2 = tracenum(2);
    try
        % Capture data from first trace
        fprintf(vsa, 'FORM:DATA REAL64');
        fprintf(vsa, sprintf(':TRAC%i:FORM "IQ"', trace1));
        flushinput(vsa);
        flushoutput(vsa);
        bufferdata1 = zeros(Ns * 2, 1);
        fprintf(vsa, sprintf('TRAC%i:DATA:XY?', trace1));
        pause(1);
        bufferdata1 = binblockread(vsa, 'double');

        % Capture data from second trace
% % % % %         fprintf(vsa, sprintf(':TRAC%i:FORM "IQ"', trace2));
% % % % %         flushinput(vsa);
% % % % %         flushoutput(vsa);
% % % % %         bufferdata2 = zeros(Ns * 2, 1);
% % % % %         fprintf(vsa, sprintf('TRAC%i:DATA:XY?', trace2));
% % % % %         pause(1);
% % % % %         bufferdata2 = binblockread(vsa, 'double');
% % % % %         
        % Convert captured data to complex format
        iqdata_vsa1 = bufferdata1(1:2:end) + 1j * bufferdata1(2:2:end);
% % % % %         iqdata_vsa2 = bufferdata2(1:2:end) + 1j * bufferdata2(2:2:end);
    catch ME
        fclose(vsa);
        rethrow(ME);
    end
end

function status_bit = check_status(bit, vsa)
    status = de2bi(mod(str2double(query(vsa, ':MEASure:STATus?')), 2^32), 32);
    status_bit = status(bit);
end