function [histBin,histValues,histBinWidth] = histparam(h)
histBin = (h.BinEdges(1:end-1) + h.BinEdges(2:end) )./2;
histValues = h.Values;
histBinWidth = h.BinEdges(2) - h.BinEdges(1);
end