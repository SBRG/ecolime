from sympy import Basic
from cobrame import mu
from collections import defaultdict


def get_reaction_stoich_dict(reaction, tol):
    stoich = {}
    for metabolite, coefficient in reaction.metabolites.items():
        if isinstance(coefficient, Basic):
            coefficient = coefficient.subs(mu, 1.)
        else:
            coefficient = coefficient

        if abs(coefficient) > tol:
            stoich[metabolite.id] = coefficient
            
    return stoich


class StoichiometryChanges(object):
    def __init__(self, old_stoich, new_stoich, tol):
        self.old_stoich, self.new_stoich = old_stoich, new_stoich
        self.set_old, self.set_new = set(old_stoich.keys()), set(
            new_stoich.keys())
        self.common = self.set_old.intersection(self.set_new)
        self.tol = tol

    @property
    def added(self):
        return self.set_new - self.common

    @property
    def removed(self):
        return self.set_old - self.common

    @property
    def changed(self):
        changed_set = set()
        for metabolite in self.common:
            diff = abs(self.old_stoich[metabolite] -
                       self.new_stoich[metabolite])
            if self.tol < diff:
                changed_set.add(metabolite)
        return changed_set


def compare_ME_reaction(old_reaction, new_reaction, output_dict, tol):
    stoich1 = get_reaction_stoich_dict(old_reaction, tol)
    stoich2 = get_reaction_stoich_dict(new_reaction, tol)

    if stoich1 != stoich2:

        changes = StoichiometryChanges(stoich1, stoich2, tol)

        if len(changes.added) > 0:
            message = 'Added Metabolites {}'.format(changes.added)
            output_dict[new_reaction.id].append(message)
        if len(changes.removed) > 0:
            message = 'Removed Metabolites {}'.format(changes.removed)
            output_dict[new_reaction.id].append(message)
        if len(changes.changed) > 0:
            for met in changes.changed:
                message = '{}: {} -> {}'.format(met,
                                                stoich1[met], stoich2[met])
                output_dict[new_reaction.id].append(message)


def find_ME_model_difference(old_model, new_model, tol):
    reaction_list = []
    output_dict = defaultdict(list)
    for old_reaction in old_model.reactions:
        try:
            new_reaction = new_model.reactions.get_by_id(old_reaction.id)
        except KeyError:
            output_dict[old_reaction.id] = ['Reaction not in new model']
            continue

        if old_reaction.lower_bound != new_reaction.lower_bound:
            output_dict[old_reaction.id] = \
                ['Reaction lower bound changed ({}->{})'.format(
                    old_reaction.lower_bound, new_reaction.lower_bound)]

        if old_reaction.upper_bound != new_reaction.upper_bound:
            output_dict[old_reaction.id] = \
                ['Reaction upper bound changed ({}->{})'.format(
                    old_reaction.upper_bound, new_reaction.upper_bound)]

        compare_ME_reaction(old_reaction, new_reaction, output_dict, tol)
        reaction_list.append(old_reaction.id)

    for new_reaction in new_model.reactions:
        if new_reaction.id not in reaction_list:
            output_dict[new_reaction.id] = ['Reaction not in old model']
    return output_dict
