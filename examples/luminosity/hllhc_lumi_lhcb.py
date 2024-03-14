from pyoptics import lumi


ip8 = lumi.IP(name="lhcb", betx=1.5, bety=1.5, px=135e-6, py=170e-6 + 1.8e-6)

b8 = lumi.Bunch(nb=2572, ips=(8))

ip8.clone(px=135e-6, py=1.8e-6 + 160e-6).luminosity(b8) / 1e38
ip8.clone(px=135e-6 - 200e-6, py=1.8e-6).luminosity(b8) / 1e38
ip8.clone(px=135e-6 + 150e-6, py=1.8e-6).luminosity(b8) / 1e38
ip8.clone(px=135e-6, py=1.8e-6 + 170e-6).luminosity(b8) / 1e38

out = []
for bet in np.linspace(0.4, 1.6, 11):
    out.append([bet, ip8.clone(betx=bet).luminosity(b8) / 1e38])


bet, l = np.array(out).T

plt.plot(bet, l, "r", label=r"$L_{\rm virt}$, $\beta^*_y=1.5$ m")
grid(True)
ylabel(r"Luminosity [$10^{34} {\rm cm^{-2}s^{-1}}$]")
xlabel(r"$\beta^*_x$ [m]")
axvline(0.5, color="k", linestyle="-")
axvline(1.5, color="k", linestyle="-")
plt.legend()
plt.tight_layout()

plt.savefig("lhcb_virt_flat.png")
