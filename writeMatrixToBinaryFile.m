function writeMatrixToBinaryFile(M,fileName)
[numRows,numCols]=size(M);
fid = fopen(fileName, 'w');
fwrite(fid,[numRows numCols], 'int32');
fwrite(fid,M,'double');
fclose(fid);
