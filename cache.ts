import { Audio } from 'expo-av';
import { StatusCodes } from 'http-status-codes';
import React, { memo, useCallback, useEffect, useMemo, useState } from 'react';
import { Dimensions, StyleSheet } from 'react-native';
import { View } from 'react-native';
import { useDispatch, useSelector } from 'react-redux';

import audioAssets, { playAudio } from '../../assets/audio/audioAssets';
import { gammaDarkGrey, gammaGrey } from '../../constants/UIConstants';
import csfConstants from '../../src/csf/csfConstants';
import dimUtils from '../../src/dimUtils';
import eyeTestResultsActions from '../../store/actions/eyeTestResultsActions';
import LandoltCCell from '../cell/LandoltCCell';
import MatrixCell from '../cell/MatrixCell';
import eyeTestUtils from '../eyeTests/eyeTestUtils';
import LoadingWheel from '../UI/LoadingWheel';

import CellGridNextButton from './CellGridNextButton';

function cosWindow2D(radius, wLength, centreX, centreY, sizeX, sizeY) {
	const pi = Math.PI;
	const cWin = Array(sizeY).fill(0).map(() => Array(sizeX).fill(0));

	for (let y = 0; y < sizeY; y++) {
		for (let x = 0; x < sizeX; x++) {
			const dist = Math.sqrt((x + 1 - centreX) ** 2 + (y + 1 - centreY) ** 2);

			if (dist <= radius - wLength / 2) {
				cWin[y][x] = 1;
			} else if (dist >= radius + wLength / 2) {
				cWin[y][x] = 0;
			} else {
				cWin[y][x] = 0.5 + 0.5 * Math.cos(pi / 2 + pi * (dist - radius) / wLength);
			}
		}
	}

	return cWin;
}

function makeFilter(filtType, fPeak, bWdth, alpha, filtSize) {
	if (typeof filtSize === 'number') {
		filtSize = [filtSize, filtSize];
	}
	const filtRadius = [Math.round(filtSize[0] / 2), Math.round(filtSize[1] / 2)];

	const madeFilter = Array(filtSize[0]).fill(0).map(() => Array(filtSize[1]).fill(0));
	const radDist = Array(filtSize[0]).fill(0).map(() => Array(filtSize[1]).fill(0));

	const X = [];
	const Y = [];
	for (let i = 0; i < filtSize[0]; i++) {
		X[i] = [];
		Y[i] = [];
		for (let j = 0; j < filtSize[1]; j++) {
			const x = j - filtRadius[1];
			const y = i - filtRadius[0];
			X[i][j] = x;
			Y[i][j] = y;
			radDist[i][j] = Math.sqrt(x * x + y * y);
		}
	}

	// Avoid divide-by-zero at the center
	radDist[filtRadius[0]][filtRadius[1]] = 0.5;

	const pi = Math.PI;

	function log2(x) {
		return Math.log(x) / Math.log(2);
	}

	for (let i = 0; i < filtSize[0]; i++) {
		for (let j = 0; j < filtSize[1]; j++) {
			const r = radDist[i][j];
			const x = X[i][j];
			const y = Y[i][j];

			if (filtType === 1) {
				madeFilter[i][j] = Math.exp(-((Math.log(2) * Math.abs(Math.log(r / fPeak)) ** 3) / ((bWdth * Math.log(2)) ** 3)));
			}
			else if (filtType === 2) {
				madeFilter[i][j] = r ** -alpha;
			}
			else if (filtType === 3) {
				let fPeakRad = fPeak * pi / 180;
				let bWdthRad = bWdth * pi / 180;
				const angDist = Math.atan2(-y, x);
				const sintheta = Math.sin(angDist);
				const costheta = Math.cos(angDist);

				const ds1 = sintheta * Math.cos(fPeakRad) - costheta * Math.sin(fPeakRad);
				const dc1 = costheta * Math.cos(fPeakRad) + sintheta * Math.sin(fPeakRad);
				let dtheta1 = Math.abs(Math.atan2(ds1, dc1));
				let value = Math.exp(((-dtheta1) ** 2) / (2 * bWdthRad ** 2));

				fPeakRad += pi;
				const ds2 = sintheta * Math.cos(fPeakRad) - costheta * Math.sin(fPeakRad);
				const dc2 = costheta * Math.cos(fPeakRad) + sintheta * Math.sin(fPeakRad);
				let dtheta2 = Math.abs(Math.atan2(ds2, dc2));
				value += Math.exp(((-dtheta2) ** 2) / (2 * bWdthRad ** 2));

				madeFilter[i][j] = value;
			}
			else if (filtType === 4) {
				const logR = log2(r);
				let val = 0.5 * (1 + Math.cos(pi * (logR - log2(fPeak))));
				if (logR > (log2(fPeak) + 1) || logR <= (log2(fPeak) - 1)) {
					val = 0;
				}
				madeFilter[i][j] = val;
			}
			else if (filtType === 5) {
				madeFilter[i][j] = Math.exp(-((log2(r) - log2(fPeak)) ** 2) / (2 * bWdth ** 2));
			}
			else if (filtType === 6) {
				madeFilter[i][j] = Math.exp(-((r - fPeak) ** 2) / (2 * bWdth ** 2));
			}
			else {
				madeFilter[i][j] = (fPeak < 0) ? (r >= Math.abs(fPeak)) : (r <= fPeak);
			}
		}
	}

	return madeFilter;
}

const CellGrid = ({
	rows,
	cols,
	screenWidthCm,
	viewDistanceCm,
	setData,
	trials,
	eyeTestResultId,
	chartId,
	setChartId,
}) => {

	const eyeTest = useSelector(state => state.eyeTests.eyeTest);

	const types = useMemo(() => eyeTestUtils.getEyeTestTypes(eyeTestResultId), [eyeTestResultId], [eyeTestResultId]);

	const dispatch = useDispatch();

	const [dims, setDims] = useState({
		container: {
			width: 0,
			height: 0
		},
		cell: {
			width: 0,
			height: 0
		},
		screenWidthCm: 0,
		viewDistanceCm,
		screenWidthDeg: 0,
		pixPerDeg: 0,
	});
	const [isGenerating, setIsGenerating] = useState(false);

	const onLayout = (event) => {
		const { scale } = Dimensions.get('window');
		const res = dimUtils.getLayout({
			event,
			cols,
			rows,
			screenWidthCm,
			viewDistanceCm,
			pixelRatio: scale,
		});
		setDims(res);
	};

	const [trialRecord, setTrialRecord] = useState();
	const [isCurrentTrialComplete, setIsCurrentTrialComplete] = useState(false);
	const [isLoading, setIsLoading] = useState(true);
	const audioPlayer = useMemo(() => new Audio.Sound(), []);

	const onNextHandler = useCallback(async () => {

		// return;
		setIsLoading(true);
		setIsCurrentTrialComplete(false);

		const updatedTrialRecord = { ...trialRecord };
		setChartId(chartId + 1);
		setIsGenerating(true);
		const strippedTrialData = { ...updatedTrialRecord };
		delete strippedTrialData.matrix;
		const nextTrialRes = await dispatch(eyeTestResultsActions.createChart({
			eocId: eyeTest.eocId,
			encounterId: eyeTest.encounterId,
			eyeTestId: eyeTest.eyeTestId,
			eyeTestResultId,
			chartId,
			data: {
				chartCount: `${trials}`,
				rowCount : `${rows.length}`,
				colCount : `${cols.length}`,
				dims,
				pixPerDeg: `${dims.pixPerDeg}`,
				trialData: strippedTrialData,
			}
		}));
		setIsGenerating(false);
		if(nextTrialRes.status !== StatusCodes.OK) {
			console.error('Error generating trial record');
			return;
		}
		if(chartId + 1 > trials ) {
			setData(eyeTestResultId);
			setTimeout(() => setIsLoading(false), 1);
			return;
		}
		const nextTrial = nextTrialRes.data.trialData;
		const couldBeRandom = nextTrialRes.data.lastPValue > csfConstants.pValueThreshold;
		const lastTrial = updatedTrialRecord.trialNo + 1 > trials;
		playAudio(audioPlayer, !lastTrial && !couldBeRandom ? audioAssets.wellDone
			: !lastTrial && couldBeRandom ? audioAssets.aimCsfWarning
				: audioAssets.oneMoreChart);
		setTrialRecord(nextTrial);
		// 1ms delay to allow the state to update before setting isLoading to false
		setTimeout(() => setIsLoading(false), 1);
	}, [trialRecord, dims, chartId, trials, dispatch, eyeTest, setData, setChartId, eyeTestResultId, rows.length, cols.length, audioPlayer]);

	useEffect(() => {
		if(isGenerating || !dims || !dims.container || !dims.container.width || !dims.container.width) return;
		const load = async () => {
			setIsGenerating(true);
			const trialRecordRes = await dispatch(eyeTestResultsActions.createChart({
				eocId: eyeTest.eocId,
				encounterId: eyeTest.encounterId,
				eyeTestId: eyeTest.eyeTestId,
				eyeTestResultId,
				chartId,
				data: {
					chartCount: `${trials}`,
					rowCount : `${rows.length}`,
					colCount : `${cols.length}`,
					dims,
					pixPerDeg: `${dims.pixPerDeg}`,
				}
			}));

			setIsGenerating(false);
			if(trialRecordRes.status !== StatusCodes.OK) {
				console.error('Error generating trial record');
				return;
			}
			const trialRecord = trialRecordRes.data.trialData;

			// const targSFcDeg = 2;
			// const pedRange = [0.001, 0.9];
			// const paramEst = [0.01, 0.5];
			// const pixPerDeg = dims.pixPerDeg;
			// const testLambda = pixPerDeg / targSFcDeg;
			// const targSizePix = dims.cellDim * 0.5;
			// const mySFFilter = makeFilter(
			// 	4,
			// 	targSizePix / testLambda,
			// 	0.5,
			// 	1,
			// 	targSizePix,
			// );
			// const mySpatialWindow = cosWindow2D(
			// 	targSizePix * 0.45,
			// 	pixPerDeg / 2,
			// 	targSizePix / 2,
			// 	targSizePix / 2,
			// 	targSizePix,
			// 	targSizePix
			// );

			// const trialRecord = {
			// 	trialNo: 0,
			// 	trialSeed: 0,
			// 	INEst: 0,
			// 	PsiEst: 0,
			// 	slopeEst: 0,
			// 	minErrEst: 0,
			// 	sigLevel: [],
			// 	targRanPedLocs: [],
			// 	pedLevel: [],
			// 	conIncLevel: [],
			// 	matchOri: [],
			// 	oriErr: [],
			// 	oriErrSorted: [],
			// 	stimSeen: [],
			// 	chartTime: -1,
			// };
			setTrialRecord(trialRecord);
			setIsLoading(false);

		};
		load();
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [dispatch, dims, rows, cols, eyeTest]);

	useEffect(() => {
		if(!trialRecord || !trialRecord.stimSeen) return;
		let isTrialComplete = true;
		for(let i = 0; i < trialRecord.stimSeen.length; i++) {
			if(trialRecord.stimSeen[i] === 0) {
				isTrialComplete = false;
				break;
			}
		}
		setIsCurrentTrialComplete(isTrialComplete);
	}, [trialRecord]);

	const styles = StyleSheet.create({
		container: {
			flex: 1,
			backgroundColor: types.isAimAcuityLowContrast ? gammaDarkGrey : gammaGrey,
			justifyContent: 'center',
			alignItems: 'center',
		},
	});

	return (
		<View
			onLayout={onLayout}
			style={styles.container}
		>
			{!isLoading && rows.map(i => {
				return (
					<View
						key={`gridCellRow${i}`}
						style={{
							flexDirection: 'row',
							justifyContent: 'center',
						}}
						testId={`gridCellRow${i}`}
					>
						{cols.map(j => {
							if(!trialRecord) return;
							const cellNumber = i * cols.length + j;

							const cellWidth = dims.cell.width / dims.pixelRatio;
							const cellHeight = dims.cell.height / dims.pixelRatio;
							const cellDim = Math.floor(Math.min(cellWidth, cellHeight));

							const matrix = trialRecord.matrix ? trialRecord.matrix[cellNumber] : [];

							return (
								<View
									style={{
										width: cellDim,
										height: cellDim,
										margin: dims.cellMargin,
										justifyContent: 'center',
										alignItems: 'center',
									}}
									key={`gridCellCol${i}x${j}`}
								>
									<View
										style={{
											width: cellDim,
											height: cellDim,
											overflow: 'hidden',
										}}
									>
										{types.isAimCsf > 0 && <MatrixCell
											trialRecord={trialRecord}
											setTrialRecord={setTrialRecord}
											cellNumber={cellNumber}
											matrix={matrix}
											cellDim={cellDim}
											eyeTestResultId={eyeTestResultId}
										/>}
										{(types.isAimAcuity || types.isAimAcuityLowContrast) && <LandoltCCell
											trialRecord={trialRecord}
											setTrialRecord={setTrialRecord}
											cellNumber={cellNumber}
											cellDim={cellDim}
											eyeTestResultId={eyeTestResultId}
										/>}

									</View>
								</View>
							);
						})}
					</View>
				);
			})}
			{isLoading && <LoadingWheel/>}
			{isCurrentTrialComplete && <CellGridNextButton onNextHandler={onNextHandler}/>}
		</View>
	);
};

export default memo(CellGrid);