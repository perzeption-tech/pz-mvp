import cv2
import cvzone
import os
from cvzone.FaceMeshModule import FaceMeshDetector
from scipy.io import savemat, loadmat
import numpy as np
import time

def measure_distance_face_cam():
    cap = cv2.VideoCapture(0)
    detector = FaceMeshDetector(maxFaces=1)
    distances = []

    # filename = 'Participant_distances_during_test.txt'
    # print(os.getcwd())
    # if os.path.exists(os.path.join(os.getcwd(), filename)):
    #     os.remove(os.path.join(os.getcwd(), filename))
    # file = open(filename, 'w')
    mat = loadmat('stop.mat')
    print(mat['stop'])
    while mat['stop'] != [[1]]:
        success, img = cap.read()
        img, faces = detector.findFaceMesh(img, draw=False)

        if faces:
            face = faces[0]
            pointLeft = face[145]
            pointRight = face[374]

            # cv2.line(img, pointLeft, pointRight, (0, 200, 0), 3)
            # cv2.circle(img, pointLeft, 5, (255, 0, 255), cv2.FILLED)
            # cv2.circle(img, pointRight, 5, (255, 0, 255), cv2.FILLED)
            w, _ = detector.findDistance(pointLeft, pointRight)
            W = 6.3

            # # Finding the Focal Length
            # d = 50
            # f = (w*d)/W
            # print(f)

            # Finding distance
            f = 515
            distance = (W * f) / w
            # print(distance)
            distances.append(distance)
            float_d = '%.2f' % distance

            cvzone.putTextRect(img, f'Distance: {float_d}cm',
                               (face[10][0] - 100, face[10][1] - 50),
                               scale=2)

        # file.write(f"{distance}\n")
        cv2.imshow("faces", img)
        savemat('test.mat', mdict={'distance': distance})
        # time.sleep(10)
        cv2.waitKey(1)
        mat = loadmat('stop.mat')
    # file.close()
