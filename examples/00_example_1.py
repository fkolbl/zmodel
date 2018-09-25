import zmodels.lti as tf
import matplotlib.pyplot as plt

H = tf.lti([1],[1,0])
print H.num
print H.den

# Test of the bode plot
w, mag, phase = tf.bode(H)

plt.figure()
plt.semilogx(w, mag)    # Bode magnitude plot
plt.grid(True, which="both")
plt.figure()
plt.semilogx(w, phase)  # Bode phase plot
plt.grid(True, which="both")

# Test of the step response
t, step = tf.step2(H)
plt.figure()
plt.plot(t, step)
plt.grid(True, which="both")

# Test of the impulse response
t, impuls = tf.impulse(H)
plt.figure()
plt.plot(t, impuls)
plt.grid(True, which="both")

## test 
H2 = tf.lti([1],[1,1])
G = H*H2
w, mag, phase = tf.bode(G)
plt.figure()
plt.semilogx(w, mag)    # Bode magnitude plot
plt.grid(True, which="both")
plt.figure()
plt.semilogx(w, phase)  # Bode phase plot
plt.grid(True, which="both")

plt.show()