function [codeword] = polar_encoder_channel_CRC(info, CRCgenerator, frozen_flag)
codeword = polar_encoder_channel(CRCgenerator(info), frozen_flag);
end

