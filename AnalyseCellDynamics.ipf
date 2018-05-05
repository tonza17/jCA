#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//////////////////////////////////////////////////////
//CALL THIS FUNCTION TO ANALYSE CALCIUM DYNAMICS  
//PARAMETERS:
//numtraces - number of cells to be analysed
//numframes - number of frames in the original movie
//minlength - only generate traces that are bigger than this value
//
//
function AnalyseCellDynamics(numtraces, numframes, minlength, framesPerSecond)
variable numtraces, numframes, minlength, framesPerSecond

wave v, tW, iW, pW,a
GetAllDynamics(v, tW, iW, pW, numTraces, minLength)
wave traces
duplicate /o traces, tr
GetAllDynamics(a, tW, iW, pW, 75, 2800)
tr[][] /= traces[p][q]
SelectTracesBiggerThanValue (tr, minlength)
wave fullSizeTRaces
RemoveIdenticalTraces (FullSizeTRaces)
killwaves traces, tr 
wave corMatrix,  results, trace, trace1, trace2, fullsizeuncortraces
killwaves corMatrix,  results, trace, trace1, trace2, fullSizeTraces, traces
SetScale/P x 0, 1/framesPerSecond,"", fullsizeuncortraces
rename fullsizeuncortraces, traces

CalculatePower (traces)
wave allPowers
allPowers[0]=nan
end

///////////////////////////////////////////////////////
function plotTrace (traceNumber)
variable traceNumber
wave traces
display /k = 1 traces[][traceNumber]


end
///////////////////////////////////////////////////////
function plotPowerSpectrum (traceNumber)
variable traceNumber
wave allPowers
display /k = 1 allPowers[][traceNumber]

end
///////////////////////////////////////////////////////
function plotWignerTransform(traceNumber, gaus, cutoffvalue)
variable traceNumber, gaus, cutoffvalue
wave traces
if (mod(DimSize(traces, 0), 2)==0)
	make /o/n=(DimSize(traces, 0)) tr
else
	make /o/n=(DimSize(traces, 0)-1) tr
endif
setscale /p x,dimoffset(traces,0),dimdelta(traces,0),waveunits(traces,0) tr
WignerTransform /Gaus=(gaus) tr
wave m_wigner
duplicate /o m_wigner, $("m_wigner"+num2str(traceNumber))
Display /k = 1;AppendImage $("m_wigner"+num2str(traceNumber))
ModifyImage $("m_wigner"+num2str(traceNumber)) ctab= {*,cutoffvalue,Rainbow,1}
killwaves m_wigner
end


////////////////////////////////////////////////////////



/////DO NOT MODIF THESE FUNCTIONS///////////////////////////////

function GetCellDynamics(w_dynamics, w_time, w_index, w_realindex, cellindex, numframes)
wave w_dynamics, w_time, w_index, w_realindex
variable cellindex, numframes
variable v_value
FindValue /v=(cellindex) w_realindex
variable index = w_index[v_Value]
duplicate /o w_index, test
//test[] = (w_index[p] == index) ? 1:0
//wavestats /q test
//v_sum = 600
make /o/n = (numframes) subset
subset = nan
variable i = 0
FindValue /V = (index) w_index
variable start = v_Value
for (i = start; i < start + numframes; i += 1)
	subset[i-start] = w_dynamics[i]
	if (w_time[i] > w_time[i + 1])
		
		break
		
	endif
endfor


end

/////////////////////////////////////////////////////////////////////////
function GetAllDynamics(w_dynamics, w_time, w_index, w_realindex, numcells, numframes)
wave w_dynamics, w_time, w_index, w_realindex
variable numcells, numframes

variable i = 0

GetCellDynamics(w_dynamics, w_time, w_index, w_realindex, 1, numframes)
wave subset

for (i = 1; i < numcells; i += 1)
	GetCellDynamics(w_dynamics, w_time, w_index, w_realindex, i, numframes)
	if (i == 1)
		duplicate /o subset, traces
	elseif (i == 2)
		concatenate {subset}, traces
	else
		concatenate /np {subset}, traces
	endif
endfor
killwaves subset
end

/////////////////////////////////////////////////////////////////////////
function SelectTracesBiggerThanValue (pop, thresholdsize)
wave pop
variable thresholdsize

variable v_mnumnans
variable numtraces = DimSize(pop, 1)
variable numpoints = DimSize(pop, 0)
make /o/n = (numpoints) trace
variable i = 0
variable found = 0
for (i = 0; i< numtraces; i +=1)
	trace[] = pop[p][i]
	wavestats /Q trace
	if (numpoints-v_numnans >= thresholdsize)
		if (found == 0)
			duplicate /o trace, fullsizetraces
		elseif (found == 1)
			concatenate {trace}, fullsizetraces
		else
			concatenate /NP {trace}, fullsizetraces
		endif
		found += 1
	endif
endfor

end

///////////////////////////////////////////////////////////////////////////
function RemoveIdenticalTraces (pop)
wave pop

variable v_mnumnans
variable numtraces = DimSize(pop, 1)
variable numpoints = DimSize(pop, 0)
make /o/n = (numpoints) trace1, trace2
variable i = 0, j = 0

//GENERATE CORRELATION MATRIX
make /o/n = (numtraces, numtraces) corMatrix
corMatrix = 0
for (i = 0; i< numtraces; i +=1)
	for (j = 0; j < numtraces; j +=1)
		trace1[] = pop[p][i]
		trace2[] = pop[p][j]
		corMatrix[i][j] = StatsCorrelation(trace1,trace2)
	endfor
endfor

//SELECT TRACES SHOWING NO CORRELATION
make /o/n =  (numtraces) results
results = 1
variable found = 0
for (i = 0; i< numtraces; i +=1)
	for (j = i; j < numtraces; j +=1)
		if ((corMatrix[i][j] ==1)&&(i != j))
			results[j] = nan
		endif
	endfor
endfor	

for (i = 0; i< numtraces; i += 1)
	trace1[] = pop[p][i]
	if (results[i] == 1)
		if (found == 0)
			duplicate /o trace1, fullsizeuncortraces
		elseif (found == 1)
			concatenate {trace1}, fullsizeuncortraces
		else
			concatenate /NP {trace1}, fullsizeuncortraces
		endif
	found += 1
	endif
endfor
end

/////////////////////////////////////////////////////////////////////

function GetSpikes(dbpop, thresholdfactor)
wave dbpop
variable thresholdfactor
variable i, j, nwaves, numpoints, deltat=1, period=3, pointnum, stddev, threshold, Pearson, totalspikes=0
wave W_FindLevels

String arrayname = "pop"
nwaves = DimSize(dbpop, 1)
numpoints= DimSize(dbpop, 0)-15
Make/O/N=(numpoints) tempwave, spikepoints
SetScale/P x 0, deltat, "", tempwave, spikepoints
Duplicate /O dbpop, spike_array, spikes
SetScale/P x 0, deltat, "", spike_array, spikes
DeletePoints /M=0 (numpoints-15), 15, spike_array, spikes
spikes=0

Make/O/N=15 W_StatsLinearCorrelationTest
Make/O/N=(nwaves) spikespercell

for(i=0; i<nwaves; i+=1)
	Duplicate/O/R=[0, numpoints-1][i] dbpop, tempwave
	//tempwave[p]=spike_array[p][i]
	Smooth/E=2 7, tempwave
	WaveStats/Q tempwave
	tempwave/=V_avg
	Differentiate/METH=2 tempwave
	WaveStats/Q tempwave
	stddev=V_sdev
	threshold=thresholdfactor*stddev
	FindLevels/P /Q/EDGE=1 tempwave, (threshold)
	Duplicate/O W_FindLevels, ycrossings
	for(j=0; j<(V_LevelsFound); j+=1)
		pointnum=(W_FindLevels[j])+1
		ycrossings[j]=tempwave[pointnum]
	endfor

	for(j=0; j<(V_LevelsFound); j+=1)
		pointnum=(W_FindLevels[j])+1
		//if((tempwave[ceil(pointnum+0.5*(period/deltat))] )<(-0.2*threshold))
			//spikes[(pointnum)][i]=0
		//elseif((tempwave[round(pointnum+1.5*(period/deltat))] )<(-0.5*threshold))
			//spikes[(pointnum)][i]=0
		//else
			spikes[(pointnum)][i]=1
			totalspikes+=V_LevelsFound
			spikespercell[i]=V_LevelsFound
		//endif
	endfor
endfor
print "totalspikes counted =" + num2str(totalspikes)
end

///////////////////////////////////////////////////////////////////////////

function countSpikes (spikepop)
wave spikepop
variable popxnum=DimSize(spikepop,0)
variable popynum=DimSize(spikepop,1)
make /o/n=(popxnum) spiketr, spiketrint
make /o/n=(popynum) spikecount
variable i=0
for (i=0; i<popynum; i+=1)
	spiketr[]=spikepop[p][i]
	Integrate spiketr/D=spiketrint
	spikecount[i]=spiketrint[popxnum]
endfor
end


function calculateSpikesDynamics(spikes, windowsize)
wave spikes
variable windowsize

variable numpoints=DimSize(spikes, 0)
variable numtr=DimSize(spikes, 1)
variable i=0, j=0, k=0

make /o/n=(numpoints) spikeDynamics
spikedynamics=0
for(i=0; i<numpoints; i+=1)
	for (j=0; j<numtr; j+=1)
		if (spikes[i][j]==1)
			spikedynamics[i]+=1
		endif
	endfor

endfor


end
///////////////////////////////////////////////////////////////////////////////
function CalculateInterspikeinterval(spikes, start, stop)
wave spikes
variable start, stop

variable numtr=DimSize(spikes, 1)
variable numpoints=DimSize(spikes, 0)

make /o/n=(numtr*numpoints) ISI
ISI=NaN

variable i=0, j=0

variable prevspike=0

for(j=0;j<numtr;j+=1)
	for(i=start; i<stop; i+=1)
		if (spikes[i][j]==1)
			ISI[i*j]=i-prevspike
			prevspike=i
			if (prevspike==stop)
				prevspike=0
			endif
		endif
	endfor
endfor

end
///////////////////////////////////////////////////////////////////////////////
function calculatePower(db)
wave db

variable numtr=DimSize(db, 1)
variable i=0
variable v_avg
variable num=0
if (mod(DimSize(db, 0), 0)==0)
	make /o/n=(DimSize(db, 0)) tr
else
	make /o/n=(DimSize(db, 0)-1) tr
endif
setscale /p x,dimoffset(db,0),dimdelta(db,0),waveunits(db,0) tr
variable start = 0
variable stop = DimSize(tr, 0)
print stop
for(i=0; i<numtr; i+=1)
	tr[]=db[p][i]
	
	duplicate /o tr, tr_FFT
	FFT/OUT=4/DEST=tr_FFT tr
	
	if (i==0)
		duplicate /o tr_fft allpowers
	elseif (i==1)
		concatenate {tr_fft}, allpowers
	else
		concatenate /NP {tr_fft}, allpowers
	endif
		num+=1
		
	
endfor
killwaves tr
end
///////////////////////////////////////////////////////////////////////////////

function spikeDynamicsPerTrace(spikes, windowsize)
wave spikes
variable windowsize

variable i = 0, j = 0 
variable numtr = DimSize (spikes, 1), numpoints = DimSize (spikes, 0)
variable v_sum

make /o/n = (numpoints/windowsize, numtr) spikefreqdynamics
make /o/n=(numpoints) trace
for (i = 0; i < numtr; i +=1)
	for (j = 0; j <  numpoints; j += windowsize)
		trace[] = spikes[p][i]
		wavestats /q/r = (j, j + windowsize) trace
		spikefreqdynamics[j/windowsize][i] = v_sum/windowsize
	endfor
endfor
end
///////////////////////////////////////////////////////////////////////////////
