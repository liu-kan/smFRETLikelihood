gamma=0.34 alph=0.08 beta=0.07

    def _calculate_fret_eff(self, pax=False):
        """Compute FRET efficiency (`E`) for each burst."""
        G = self.get_gamma_array()
        if not pax:
            E = [na / (g * nd + na) for nd, na, g in zip(self.nd, self.na, G)]
        else:
            adr = self._aex_dex_ratio
            E = [na * (1 + adr) / (g * (nd + nda) + na * (1 + adr))
                 for nd, na, nda, g in zip(self.nd, self.na, self.nda, G)]
        self.add(E=E, pax=pax)

    def _calculate_stoich(self, pax=False):
        """Compute "stoichiometry" (the `S` parameter) for each burst."""
        G = self.get_gamma_array()
        naa = self.naa
        if 'PAX' in self.meas_type:
            naa = self.naa_aexonly
        if not pax:
            S = [(g * nd + na) / (g * nd + na + naa / self.beta)
                 for nd, na, naa, g in zip(self.nd, self.na, naa, G)]
        else:
            # This is a PAX-enhanced formula which uses information
            # from both alternation periods in order to compute S
            alpha = 1 - self._aex_fraction
            S = [(g * (nd + nda) + na / alpha) /
                 (g * (nd + nda) + na / alpha + naa / (alpha * self.beta))
                 for nd, na, nda, naa, g in
                 zip(self.nd, self.na, self.nda, naa, G)]
        self.add(S=S)

Eco=(x*nb-ba)/(nb-ba-bd)

    nb=5938.650324780
    ba=995.0
    bd=3237.39

