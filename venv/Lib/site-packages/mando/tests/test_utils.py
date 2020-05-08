import pytest
from mando.utils import action_by_type, ensure_dashes, find_param_docs, split_doc


ACTION_BY_TYPE_CASES = [
    (True, {'action': 'store_false'}),
    (False, {'action': 'store_true'}),
    ([], {'action': 'append'}),
    ([1, False], {'action': 'append'}),
    (None, {}),
    (1, {'type': type(1)}),
    (1.1, {'type': type(1.1)}),
    ('1', {'type': type('1')}),
]


@pytest.mark.parametrize('obj,result', ACTION_BY_TYPE_CASES)
def test_action_by_type(obj, result):
        assert result == action_by_type(obj)


ENSURE_DASHES_CASES = [
    (['m'], ['-m']),
    (['m', 'min'], ['-m', '--min']),
    (['-m'], ['-m']),
    (['-m', 'min'], ['-m', '--min']),
    (['m', '--min'], ['-m', '--min']),
    (['-m', '--min'], ['-m', '--min']),
    (['-m', '--min', 'l', 'less'], ['-m', '--min', '-l', '--less']),
]


@pytest.mark.parametrize('opts,result', ENSURE_DASHES_CASES)
def test_ensure_dashes(opts, result):
    assert result == list(ensure_dashes(opts))


SPLIT_DOC_CASES = [
    ('', ['', '']),
    ('only help.', ['only help.', 'only help.']),
    ('help.\nstill help.', ['help.\nstill help.', 'help.\nstill help.']),
    ('help\n\ndesc', ['help', 'desc']),
    ('help\n\n\ndesc\n', ['help', 'desc']),
]


@pytest.mark.parametrize('doc,parts', SPLIT_DOC_CASES)
def test_split_doc(doc, parts):
    assert parts == split_doc(doc)


a_1 = {'a_param': (['a-param'], {'help': 'Short story.'})}
a_1_1 = {'a_param': (['a_param'], {'help': 'Short story.'})}
a_2 = {'j': (['-j'], {'help': 'Woow'})}
a_3 = {'noun': (['-n', '--noun'], {'help': 'cat'})}
a_all = {}
for a in (a_1, a_2, a_3):
    a_all.update(a)


FIND_PARAM_CASES = [
    ('', {}),
    ('Brevity is the soul of wit.', {}),
    (':param a-param: Short story.', a_1),
    (':param a_param: Short story.', a_1_1),
    (':param -j: Woow', a_2),
    (':param -n, --noun: cat', a_3),
    ('''
         Some short text here and there.

         :param well: water''', {'well': (['well'], {'help': 'water'})}),
    ('''
         :param a-param: Short story.
         :param -j: Woow
         :param -n, --noun: cat''', a_all),
    ('''
         Lemme see.

         :param long-story: A long story believe me: when all started, Adam and Bob were just two little farmers.
         ''', {'long_story': (['long-story'], {'help': 'A long story '
                              'believe me: when all started, Adam and '
                              'Bob were just two little farmers.'})}),
]


@pytest.mark.parametrize('doc,params', FIND_PARAM_CASES)
def test_find_param(doc, params):
        found_params = find_param_docs(doc)
        assert params.keys() == found_params.keys()
        for key, value in params.items():
            assert key in found_params
            found_value = found_params[key]
            assert value[0] == found_value[0]
            for kwarg, val in value[1].items():
                assert val == found_value[1][kwarg]
