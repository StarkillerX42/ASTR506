program solsys
    use integrators
    implicit none
    real, dimension(0:9) :: masses
    real, dimension(0:9) :: sma
    real, dimension(0:9) :: ecc
    real, dimension(0:9) :: incl
    real, dimension(0:9) :: bigomega
    real, dimension(0:9) :: omega
    real, dimension(0:9) :: ftrue

 ! masses in g, sma in AU, angles in deg
 ! Orbital parameters and masses from dePater & Lissauer for planets
 ! Pluto info from
 ! https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
 masses = (/ 1.989e33, 0.3302e27, 4.8685e27, 5.9736e27, 0.64185e27, &
             1898.6e27, 568.46e27, 86.832e27, 102.43e27, 0.01303e27 /)
 sma = (/ 0.0, 0.38709880, 0.72333201, 1.00000083, 1.52368946, &
          5.2027584, 9.5428244, 19.19206, 30.06893, 39.48168677 /)
 ecc = (/ 0.0, 0.20563175, 0.00677177, 0.016708617, 0.09340062, &
          0.048495, 0.055509, 0.04630, 0.00899, 0.24880766 /)
 incl = (/ 0.0, 7.00499, 3.39447, 0.0, 1.84973, 1.3033, 2.4889, &
           0.773, 1.770, 17.14175 /)
 bigomega = (/ 0.0, 48.3309, 76.6799, 0.0, 49.5581, 100.464, &
               113.666, 74.01, 131.787, 110.30347 /)
 omega = (/ 0.0, 77.4561, 131.5637, 102.9374, 336.6023, &
            14.331, 93.057, 173.01, 48.12, 224.06676 /)
 
 
 print *, masses(2), incl(5)
end program solsys                     
