"""Module for writing Chemkin-format mechanism."""

from textwrap import fill

from .chem_utilities import units, GAS_CONSTANT, get_elem_wt
from .mech_interpret import elem_wt


def write_mech(filename, elems, specs, reacs):
    """Write Chemkin-format mechanism.

    Input
    """
    with open(filename, 'w') as file:

        # elements
        file.write('ELEMENTS\n')

        elem_wt_orig = get_elem_wt()
        elem_new = set(elem_wt.items()) - set(elem_wt_orig.items())
        elem_new = dict(elem_new)
        for e in elems:
            # write atomic weight if necessary
            if e in elem_new:
                file.write(e + ' /' + str(elem_wt[e.lower()]) + '/ \n')
            else:
                file.write(e + '\n')

        file.write('END\n\n')

        # write list of species
        species_names = fill(
            '  '.join([sp.name for sp in specs]),
            width=60,
            break_long_words=False,
            break_on_hyphens=False
            )
        file.write(
            'SPECIES\n' + 
            f'{species_names}\n'
            'END\n\n'
            )

        # reactions
        file.write('REACTIONS                   CAL/MOLE\n')

        for rxn in reacs:
            line = ''

            # reactants
            for sp in rxn.reac:
                isp = rxn.reac.index(sp)
                # print stoich coefficient if other than one
                if rxn.reac_nu[isp] != 1:
                    line += str(rxn.reac_nu[isp]) + sp
                else:
                    line += sp

                if (len(rxn.reac) - 1) > isp:
                    line += '+'

            # third body in reactants
            if rxn.pdep:
                if rxn.pdep_sp:
                    line += '(+{:s})'.format(rxn.pdep_sp)
                else:
                    line += '(+m)'
            elif rxn.thd_body:
                line += '+m'

            if rxn.rev:
                line += '='
            else:
                line += '=>'

            # products
            for sp in rxn.prod:
                isp = rxn.prod.index(sp)
                # print stoich coefficient if other than one
                if rxn.prod_nu[isp] != 1:
                    line += str(rxn.prod_nu[isp]) + sp
                else:
                    line += sp

                if (len(rxn.prod) - 1) > isp:
                    line += '+'

            # third body in products
            if rxn.pdep:
                if rxn.pdep_sp:
                    line += '(+{:s})'.format(rxn.pdep_sp)
                else:
                    line += '(+m)'
            elif rxn.thd_body:
                line += '+m'

            # Convert internal units to moles
            reac_ord = sum(rxn.reac_nu)
            pre_factor = rxn.rate_parameters.pre_factor
            if rxn.thd_body:
                pre_factor *= 1000. ** reac_ord
            elif rxn.pdep:
                # Low- (chemically activated bimolecular reaction) or
                # high-pressure (fall-off reaction) limit parameters
                pre_factor *= 1000. ** (reac_ord - 1.)
            else:
                # Elementary reaction
                pre_factor *= 1000. ** (reac_ord - 1.)

            # now add Arrhenius coefficients to the same line
            line += ' {:.4e} {:.4e} {:.4e}'.format(
                pre_factor, 
                rxn.rate_parameters.temp_exponent, 
                (rxn.rate_parameters.act_energy * GAS_CONSTANT).to(units('cal/mole')).magnitude
                )

            line += '\n'
            file.write(line)

            # line for reverse Arrhenius parameters, if any
            if rxn.rev and rxn.rev_par:
                # Convert internal units to moles
                reac_ord = sum(rxn.prod_nu)
                pre_factor = rxn.rev_par.pre_factor
                if rxn.thd_body:
                    pre_factor *= 1000. ** reac_ord
                elif rxn.pdep:
                    # Low- (chemically activated bimolecular reaction) or
                    # high-pressure (fall-off reaction) limit parameters
                    pre_factor *= 1000. ** (reac_ord - 1.)
                else:
                    # Elementary reaction
                    pre_factor *= 1000. ** (reac_ord - 1.)

                line = '  rev/ {:.4e}  {:.4e}  {:.4e} /\n'.format(
                    pre_factor,
                    rxn.rev_par.temp_exponent,
                    (rxn.rev_par.act_energy * GAS_CONSTANT).to(units('cal/mole')).magnitude
                    )
                file.write(line)

            # write Lindemann low- or high-pressure limit Arrhenius parameters
            if rxn.pdep:
                if rxn.low:
                    line = '  low /{:.4e}  {:.4e}  {:.4e} /\n'.format(
                        rxn.low.pre_factor * 1000. ** sum(rxn.reac_nu), 
                        rxn.low.temp_exponent, 
                        (rxn.low.act_energy * GAS_CONSTANT).to(units('cal/mole')).magnitude
                        )
                else:
                    line = '  high /{:.4e}  {:.4e}  {:.4e} /\n'.format(
                        rxn.high.pre_factor * 1000. ** (sum(rxn.reac_nu) - 2.), 
                        rxn.high.temp_exponent, 
                        (rxn.high.act_energy * GAS_CONSTANT).to(units('cal/mole')).magnitude
                        )
                file.write(line)

            # write Troe parameters if any
            if rxn.troe:
                troe = rxn.troe_par
                if len(troe) == 3:
                    line = '  troe/ {:.4e} {:.4e} {:.4e} /\n'.format(
                        troe[0], troe[1], troe[2]
                        )
                else:
                    line = '  troe/ {:.4e} {:.4e} {:.4e} {:.4e} /\n'.format(
                        troe[0], troe[1], troe[2], troe[3]
                        )
                file.write(line)

            # write SRI parameters if any
            if rxn.sri:
                sri = rxn.sri_par
                if len(sri) == 3:
                    line = '  sri/ {:.4e} {:.4e} {:.4e} /\n'.format(
                        sri[0], sri[1], sri[2]
                        )
                else:
                    line = '  sri/ {:.4e} {:.4e} {:.4e} {:.4e} {:.4e} /\n'.format(
                        sri[0], sri[1], sri[2], sri[3], sri[4]
                        )
                file.write(line)

            # write CHEB parameters, if any
            if rxn.cheb:
                line = (
                    '  pcheb / {:.2f} '.format(rxn.cheb_plim[0].to('atm').magnitude) +
                    '{:.2f} /\n'.format(rxn.cheb_plim[1].to('atm').magnitude) +
                    '  tcheb / {:.1f} '.format(rxn.cheb_tlim[0]) +
                    '{:.1f} /\n'.format(rxn.cheb_tlim[1]) +
                    '  cheb / {} {} '.format(rxn.cheb_n_temp, rxn.cheb_n_pres)
                    )
                file.write(line)

                line = '  cheb /'
                for par in rxn.cheb_par:
                    if len(line) > 70:
                        file.write(line + ' /\n')
                        line = '  cheb /'
                    line += ' {: 7.5e}'.format(par)
                file.write(line + line + ' /\n')

            # write PLOG parameters, if any
            if rxn.plog:
                for par in rxn.plog_par:
                    line = (
                        '  plog/ {:.2e} '.format(par[0].to('atm').magnitude) +
                        '{:.4e} {:.4e} {:.4e} /\n'.format(
                            par[1].pre_factor * 1000. ** (sum(rxn.reac_nu) - 1.),
                            par[1].temp_exponent,
                            (par[1].act_energy * GAS_CONSTANT).to(units('cal/mole')).magnitude
                            )
                        )
                    file.write(line)

            # third-body efficiencies
            if len(rxn.thd_body_eff) > 0:
                line = '  '
                for thd_body in rxn.thd_body_eff:
                    thd_eff = '{:.2f}'.format(thd_body[1])
                    line += thd_body[0] + '/' + thd_eff + '/ '

                    # move to next line if long
                    if (len(line) >= 60 and
                        (rxn.thd_body_eff.index(thd_body)
                        is not (len(rxn.thd_body_eff)-1)
                        )
                        ):
                        line += '\n'
                        file.write(line)
                        line = '  '

                line += '\n'
                file.write(line)

            # duplicate reaction flag
            if rxn.dup:
                file.write('  DUPLICATE\n')

        file.write('END')
