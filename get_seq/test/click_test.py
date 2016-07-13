import click
from click.testing import CliRunner
# py.test click.test.py

def testing_test():
    @click.command()
    @click.argument('x',nargs=1)
    @click.option('--plus_two/--no', default='n')
    def test(x, plus_two):
        x = int(x)
        if not plus_two:
            click.echo(x)
        if plus_two:
            click.echo(x+2)


    runner = CliRunner()
    result = runner.invoke(test, ["2"])
    result2 = runner.invoke(test, ["2", "--plus_two"])
    assert result.exit_code == 0
    assert result.output == '2\n'
    assert result2.exit_code == 0
    assert result2.output == '4\n'

if __name__ == '__main__':
    testing_test()
